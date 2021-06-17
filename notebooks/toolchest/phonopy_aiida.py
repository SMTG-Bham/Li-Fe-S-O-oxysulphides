"""
Tool for working with Phonopy calculations
"""
from pathlib import Path
from warnings import warn
import math
from itertools import product
import numpy as np
from pymatgen.core import Structure

from aiida.orm import StructureData, Dict
import aiida.orm as orm
from aiida.engine import calcfunction
from ase import Atoms

from aiida_vasp.parsers.file_parsers.poscar import PoscarParser
from aiida.common.exceptions import NotExistent

from sumo.plotting import phonon_bs_plotter
import yaml


def get_phonon_obj(work, nac='auto', settings_update=None):
    """
    Reconstruct a phonopy object from finished/unfinished workchain.

    Useful when recomputing certain properties are needed.
    """
    force_set = work.outputs.force_sets
    if nac == 'auto':
        if 'nac_params' in work.outputs:
            params = {'nac_params': work.outputs.nac_params}
    elif nac:
        params = {'nac_params': nac}
    else:
        params = {}

    settings = work.outputs.phonon_setting_info.get_dict()

    if settings_update:
        settings.update(settings_update)
    if 'relaxed_structure' in work.outputs:
        phonon = get_phonopy_instance(work.outputs.relaxed_structure, settings, params)
    else:
        # No relaxation - use the input structure directly
        phonon = get_phonopy_instance(work.inputs.structure, settings, params)

    # Treat the magmom
    try:
        phonon.unitcell.set_magnetic_moments(work.inputs.structure.get_attribute('MAGMOM'))
    except (KeyError, AttributeError):
        pass

    displacements = work.outputs.phonon_setting_info['displacement_dataset']
    phonon.dataset = displacements
    phonon.set_forces(force_set.get_array('force_sets'))
    try:
        phonon.set_force_constants(work.outputs.force_constants.get_array('force_constants'))
    except (AttributeError, KeyError, NotExistent):
        phonon.produce_force_constants()
        warn('Cannot locate force constants - producing force constants from force_sets.')
    return phonon

def get_phonopy_instance(structure, phonon_settings_dict, params):
    from phonopy import Phonopy
    from aiida_phonopy.common.utils import phonopy_atoms_from_structure

    unit_cell = phonopy_atoms_from_structure(structure)

    # Convert to plain python dictionary
    if isinstance(phonon_settings_dict, Dict):
        phonon_settings_dict = phonon_settings_dict.get_dict()

    # Set the magnetic moments
    has_magmom = 'magmom' in phonon_settings_dict
    if has_magmom:
        unit_cell.set_magnetic_moments(phonon_settings_dict['magmom'])
        primitive_matrix = phonon_settings_dict.get('primitive_matrix', None)
    else:
        primitive_matrix = phonon_settings_dict.get('primitive_matrix', 'auto')
        
    phonon = Phonopy(
        unit_cell,
        supercell_matrix=phonon_settings_dict['supercell_matrix'],
        primitive_matrix=primitive_matrix,
        symprec=phonon_settings_dict['symmetry_tolerance'])
    if 'nac_params' in params:
        from phonopy.interface.calculator import get_default_physical_units
        units = get_default_physical_units('vasp')
        factor = units['nac_factor']
        nac_params = {'born': params['nac_params'].get_array('born_charges'),
                      'dielectric': params['nac_params'].get_array('epsilon'),
                      'factor': factor}
        phonon.nac_params = nac_params

    return phonon

def mode_mapping_gamma_from_work_node(work, qstart, qfinish, qsample, band, dryrun=True):
    """
    Generate mode mapping using a work node
    """
    from aiida.orm import Node, Int, Float
    force_set = work.outputs.force_sets
    if not isinstance(band, Node):
        band = Int(band)
    if not isinstance(qstart, Node):
        qstart = Float(qstart)
    if not isinstance(qfinish, Node):
        qfinish = Float(qfinish)
    if not isinstance(qsample, Node):
        qsample = Int(qsample)

    if 'relaxed_structure' in work.outputs:
        args = [
            work.outputs.relaxed_structure, work.outputs.phonon_setting_info, force_set, work.outputs.nac_params, qstart, qfinish, qsample,
            band
        ]
    else:
        # No relaxation - use the input structure directly
        args = [work.inputs.structure, work.outputs.phonon_setting_info, force_set, work.outputs.nac_params, qstart, qfinish, qsample, band]
    if dryrun:
        kwargs = {'metadata': {'store_provenance': False}}
    else:
        kwargs = {}

    return mode_mapping_gamma(*args, **kwargs)


@calcfunction
def mode_mapping_gamma(structure, phonon_settings, force_set, nac_param, qstart, qfinish, qsamples, band):
    """Generate pushed structures at gamma point"""
    phonon = get_phonopy_instance(
        structure,
        phonon_settings,
        {'nac_params': nac_param},
    )
    displacements = phonon_settings['displacement_dataset']
    phonon.dataset = displacements
    phonon.set_forces(force_set.get_array('force_sets'))
    phonon.produce_force_constants()

    frames = {}
    qscale = math.sqrt(len(phonon.unitcell.get_scaled_positions()))
    qpoints = []
    for i, q in enumerate(np.linspace(qstart.value, qfinish.value, qsamples.value)):
        phonon.set_modulations((1, 1, 1), [[[0, 0, 0], band.value, q * qscale, 0.0]])
        cell = phonon.get_modulated_supercells()[0]
        atoms = Atoms(positions=cell.positions, cell=cell.cell, numbers=cell.numbers, pbc=True)
        struct = StructureData(ase=atoms)
        struct.label = structure.label + f' q_{i:03d}'
        frames[f'q_{i:03d}'] = struct
        qpoints.append(q)
    frames['info'] = Dict(dict={'Q_list': qpoints, 'band': band.value, 'qscale': qscale})
    return frames


@calcfunction
def mode_mapping_1d(structure: orm.StructureData, phonon_settings: orm.Dict, force_set: orm.ArrayData, nac_param: orm.Dict,
                    qstart: orm.Float, qfinish: orm.Float, qsamples: orm.Int, band: orm.Int, q_point: orm.List, supercell: orm.List):
    """
    Generate pushed structures at any qpoint

    Args:
        structure: the input structure
        phonon_settings: the settings of the phonopy
        force_set: computed force_set for phonopy
        nac_param: NAC parameters used for Phonopy
        qstart: The start push amplitude
        qfinish: The finish push amplitude
        qsample: Number of samples for push
        band: Id of the band at q_point to be pushed
        q_point: The qpoint at which the mode should be pushed
        sueprcell: The supercell expansion for which the mode map to be calculated.

    Returns:
        A dictionary of pushed frames and mode mapping information.

    """
    phonon = get_phonopy_instance(
        structure,
        phonon_settings,
        {'nac_params': nac_param},
    )
    displacements = phonon_settings['displacement_dataset']
    phonon.dataset = displacements
    phonon.set_forces(force_set.get_array('force_sets'))
    phonon.produce_force_constants()

    frames = {}
    qscale = math.sqrt(len(phonon.unitcell.get_scaled_positions()))
    qpoints = []
    for i, q in enumerate(np.linspace(qstart.value, qfinish.value, qsamples.value)):
        phonon.set_modulations(supercell.get_list(), [[q_point.get_list(), band.value, q * qscale, 0.0]])
        cell = phonon.get_modulated_supercells()[0]
        atoms = Atoms(positions=cell.positions, cell=cell.cell, numbers=cell.numbers, pbc=True)
        struct = StructureData(ase=atoms)
        struct.label = structure.label + f' q_{i:03d}'
        frames[f'q_{i:03d}'] = struct
        qpoints.append(q)
    frames['info'] = Dict(dict={'Q_list': qpoints, 'band': band.value, 'qscale': qscale})
    return frames


@calcfunction
def mode_mapping_gamma_2d(structure, phonon_settings, force_set, nac_param, settings):
    """
    Generate pushed structures at gamma point

    :param: settings - a dictionary with "qlist1", "qlist2" and "band1", "band2"
    """
    phonon = get_phonopy_instance(
        structure,
        phonon_settings,
        {'nac_params': nac_param},
    )
    displacements = phonon_settings['displacement_dataset']
    phonon.dataset = displacements
    phonon.set_forces(force_set.get_array('force_sets'))
    phonon.produce_force_constants()

    frames = {}
    qscale = math.sqrt(len(phonon.unitcell.get_scaled_positions()))
    qpoints = []
    qlist1 = settings['qlist1']
    qlist2 = settings['qlist2']
    band1, band2 = settings['band1'], settings['band2']
    for i, (q1, q2) in enumerate(product(qlist1, qlist2)):
        phonon.set_modulations((1, 1, 1), [[[0, 0, 0], band1, q1 * qscale, 0.0], [[0, 0, 0], band2, q2 * qscale, 0.0]])
        phonon.write_modulations()
        struct = StructureData(pymatgen=Structure.from_file('MPOSCAR'))
        struct.label = structure.label + f' disp_{i:03d}'
        frames[f'disp_{i:03d}'] = struct
        qpoints.append([q1, q2])
    frames['info'] = Dict(dict={'Q_list': qpoints, 'band1': band1, 'band2': band2, 'qscale': qscale})
    return frames


### SUMO interface ####

class SumoPlotInterface:
    """Plot phonon objects from Phonopy using sumo"""
    default_plot_args_phonon = {
        'units':"Thz",
        'height': 6.0,
        'width': 6.0,
        'style': None,
        'no_base_style': False,
        'from_json': None,
        'legend': None,
    }
    def __init__(self, structure: orm.StructureData, phonon_obj=None, symprec=1e-5):
        """
        Plotting interface based on phonon object using sumo
        """
        self.phonon_obj = phonon_obj
        self.structure = structure

        # initialise the attributes
        self.kpath = None
        self.kpoints = None
        self.labels = None
        self.bs = None
        self.dos = None
        self.symprec = symprec

        # Compute the band path
        self.band_path = self.get_band_path()

    def get_band_path(self, mode='bradcrack', symprec=None, spg=None, line_density=60, cart_coords=False, kpt_list=None, labels=None, phonopy=None):
        """
        Compute the band pathway using sumo, wraps around the sumo.symmetry.kpoints.get_band_path function
        """
        from sumo.symmetry.kpoints import get_path_data
        self.kpath, self.kpoints, self.labels = get_path_data(
        self.structure.get_pymatgen(),
        mode=mode,
        symprec=self.symprec if symprec is None else symprec,
        spg=spg,
        kpt_list=kpt_list,
        labels=labels,
        phonopy=True,
        line_density=line_density,
    )
        return self.kpath, self.kpoints, self.labels

    def get_phonon_band_structure(self, include_eigenvectors=False, savepath=None):
        """Compute the actual band structure using phonopy based on the band structure comptued"""
        import tempfile
        import shutil
        from pymatgen.io.phonopy import get_ph_bs_symm_line
        self.phonon_obj.set_band_structure(self.kpoints, is_eigenvectors=include_eigenvectors, labels=self.labels)
        if savepath is None:
            tempd = tempfile.mkdtemp()
            yaml_file = str(Path(tempd) / 'sumo_bands.yaml')
        else:
            yaml_file = savepath
        self.phonon_obj._band_structure.write_yaml(filename=yaml_file)
        self.bs = get_ph_bs_symm_line(yaml_file, has_nac=False, labels_dict=self.kpath.kpoints)
        # Remove the temporary file
        if savepath is None:
            shutil.rmtree(tempd)
        return self.bs

    def get_dos(self, qmesh=(8, 8, 8)):
        """
        Compute the dos from phonopy
        """
        self.phonon_obj.set_mesh(
            qmesh, 
            is_gamma_center=False,
            is_eigenvectors=True,
            is_mesh_symmetry=False,
        )
        self.phonon_obj.set_total_DOS()
        dos_freq, dos_val = self.phonon_obj.get_total_DOS()
        dos = np.zeros((len(dos_freq), 2))
        dos[:, 0], dos[:, 1] = dos_freq, dos_val
        self.dos = dos
        return dos

    def get_plotter(self):
        from sumo.plotting.phonon_bs_plotter import SPhononBSPlotter
        plotter = SPhononBSPlotter(self.bs)
        return plotter
        
    def plot_phonon_bs(self, plotter=None, save_prefix=None, image_format='pdf', dpi=200, **kwargs):
        """
        Plot the phonon band structure, note that the phonon band structure must be computed
        first
        """
        from matplotlib import rcParams
        from sumo.cli.phonon_bandplot import save_data_files

        if plotter is None:
            plotter = self.get_plotter()

        extra_args = dict(self.default_plot_args_phonon)
        extra_args.update(**kwargs)
        plt = plotter.get_plot(
            dos=self.dos,
            **extra_args
        )

        if save_prefix:
            save_prefix = Path(save_prefix)
            directory = str(save_prefix.parent)
            prefix = save_prefix.name
            basename = "{}_phonon_band.{}".format(save_prefix, image_format)
            filename = "{}_{}".format(prefix, basename) if prefix else basename

            if dpi is None:
                dpi = rcParams["figure.dpi"]
            plt.savefig(filename, format=image_format, dpi=dpi, bbox_inches="tight")

            filename = save_data_files(plotter.bs, prefix=prefix, directory=directory)
            return filename
        else:
            return plt


    @classmethod
    def from_work_node(cls, work, settings_udpate=None, **kwargs):
        if 'relaxed_structure' in work.outputs:
            structure = work.outputs.relaxed_structure
        else:
            structure = work.inputs.structure

        return cls(structure, phonon_obj=get_phonon_obj(work, settings_update=settings_udpate), **kwargs)
