"""
Tools for assisting convex hull construction using AiiDA/atomate data
"""
from tqdm.auto import tqdm
from typing import List, Union

import pandas as pd

from aiida_vasp.parsers.file_parsers.potcar import MultiPotcarIo
from pymatgen.entries.computed_entries import ComputedEntry
from pymatgen.core import Structure

from pymatgen.entries.compatibility import Correction
from collections import defaultdict

# Import all AiiDA models
import aiida.orm as orm
from aiida.orm import Group, StructureData, KpointsData, QueryBuilder, Dict, Str, Float, Node
from aiida.plugins import WorkflowFactory, CalculationFactory

_pmg_s_cache = {}


def print_spg_with_limit(pmg_struct, tol_start=-1, tol_finish=-9):
    """Print the Space Group symbol with corresponding tolerance"""
    for tol in range(tol_start, tol_finish, -1):
        spg, _ = pmg_struct.get_space_group_info(10 ** (tol), -1)
        print(f'{spg:<10} 1.0e{tol}')

def get_pmg_struct(node: orm.StructureData) -> Structure:
    """
    Caching for pymatgen structures
    """
    if node.uuid in _pmg_s_cache:
        return _pmg_s_cache[node.uuid]
    else:
        ps = node.get_pymatgen()
        _pmg_s_cache[node.uuid] = ps
        return ps


def get_u_elem(struc, ldauu, elem):
    """
    Reliabily get the value of U for Fe.
    Return -1 if the entry does not have Fe - so compatible with any U calculations
    """
    species = MultiPotcarIo.potentials_order(struc)
    if elem in species:
        ife = species.index(elem)
        if ldauu is None:
            return 0.
        return ldauu[ife]
    return -1


def get_u_map(struc: orm.StructureData, ldauu: List[int]) -> dict:
    """
    Reliabily get the value of U for Fe.
    Return -1 if the entry does not have Fe - so compatible with any U calculations
    """
    species = MultiPotcarIo.potentials_order(struc)
    mapping = {}
    for symbol in species:
        isym = species.index(symbol)
        if ldauu is None:
            mapping[symbol] = 0.0
        else:
            mapping[symbol] = ldauu[isym]
    return mapping


def filter_umap(umap_series: pd.Series, target_map: dict) -> Union[bool]:
    """
    Return a mask of valid records for a specific U assignment
    """
    mask = []
    for umap in umap_series:
        ok = True
        for elem, value in target_map.items():
            if elem in umap and umap[elem] != value:
                ok = False
        mask.append(ok)
    return mask


def get_functional(incar: dict, pot: str) -> str:
    """
    Return the name of the functional

    Args:
        incar (dict): A dictionary for setting the INCAR
        pot (str): Potential family
    """
    if incar.get('metagga'):
        return incar.get('metagga').lower()

    if pot.startswith('LDA'):
        if incar.get('gga'):
            return 'gga+ldapp'
        else:
            return 'lda'
    elif pot.startswith('PBE'):
        gga = incar.get('gga')
        hf = incar.get('lhfcalc')
        if not hf:
            if (not gga) or gga.lower() == 'pe':
                return 'pbe'
            if gga.lower() == 'ps':
                return 'pbesol'
        else:
            if (not gga) or gga.lower() == 'pe':
                if incar.get('aexx') in [
                        0.25, None
                ] and (incar.get('hfscreen') - 0.2 < 0.01):
                    return 'hse06'

    return 'unknown'


def get_entry(df: pd.DataFrame,
              pmg_col='pmg_struct',
              label_col='label',
              uuid_col='structure_uuid',
              umap_col='umap',
              xc_col='functional',
              eng_col='energy') -> List[ComputedEntry]:
    """
    Create entry from the dataframe containing data, typically returned from the
    `get_relax_records` function.
    """
    pd_entries = []
    for idx, row in df.iterrows():
        comp = row[pmg_col].composition
        attrs = {
            'struct_name': row[label_col],
            'entry_type': 'MP' if 'mp' in row[label_col] else 'AIRSS',
            'structure_uuid': row[uuid_col],
            'calc_u': row[umap_col],
            'functional': row[xc_col],
            'volume': row[pmg_col].volume, 
            'dataframe_idx': idx,  # Store the entry index in the original dataframe
        }
        pd_entries.append(
            ComputedEntry(comp, energy=row[eng_col], parameters=attrs))
    return pd_entries


def get_relax_record_single(node: orm.WorkChainNode) -> pd.DataFrame:
    """
    Obtained the relaxation record for a single node - this returns a DataFrame that
    can be used as the input for the `get_entry` function.
    """
    q = QueryBuilder()
    q.append(Node,
             filters={'id': node.id},
             tag='root',
             project=['label', 'uuid'])  # The relaxation workchain
    q.append(StructureData, with_outgoing='root',
             project=['*', 'label'])  # The structure
    q.append(Str,
             with_outgoing='root',
             edge_filters={'label': {
                 'like': '%potential_family'
             }},
             project=['attributes.value'])  # PP
    q.append(
        Dict,
        with_outgoing='root',
        edge_filters={'label': {
            'and': [
            {'like': '%parameters'},
            {'!like': 'static%'},
            ]
        }},
        project=['attributes'],
    )
    q.append(Dict,
             with_incoming='root',
             edge_filters={'label': 'misc'},
             project=['attributes.total_energies', 'attributes.maximum_force'
                      ])  # support for both new and old style
    q.append(StructureData, with_incoming='root',
             project=['*'])  # output structure for the relax
    records = []
    for relax_label, relax_uuid, struct, label, pp, inpd, eng, fmax, relaxed_structure in q.iterall(
    ):
        # Deconstruct as some of the structures does not have vasp in the input
        if 'vasp' in inpd:
            inpd = inpd['vasp']
        else:
            inpd = inpd['incar']

        # Unpack the energy
        if 'energy_no_entropy' in eng:
            eng = eng['energy_no_entropy']
        else:
            eng = eng['energy_extrapolated']
        xc = get_functional(inpd, pp)

        pstruct_in = get_pmg_struct(struct)  # This is the input structure
        pstruct_out = get_pmg_struct(relaxed_structure)
        spg_in, _ = pstruct_in.get_space_group_info(1e-5, -1)
        spg_out, _ = pstruct_out.get_space_group_info(1e-5, -1)
        spg_in_high_tol, _ = pstruct_in.get_space_group_info(1e-1, -1)
        spg_out_high_tol, _ = pstruct_out.get_space_group_info(1e-1, -1)

        records.append({
            'group_label':
            None,
            'original_label':
            label,
            'original_uuid':
            struct.uuid,
            'structure':
            struct,
            'structure_uuid':
            struct.uuid,
            'label':
            label,
            'relax_label':
            relax_label,
            'relax_uuid':
            relax_uuid,
            'energy':
            eng,
            'functional':
            xc,
            'pp':
            pp,
            'max_force':
            fmax,
            'umap':
            get_u_map(struct, inpd.get('ldauu')),
            'note':
            'relax',
            'relax_structure':
            relaxed_structure,
            'pmg_struct':
            pstruct_in,
            'formula':
            pstruct_in.composition.reduced_formula,
            'volume':
            pstruct_out.volume,
            'volume':
            pstruct_out.volume,
            'volume_per_fu':
            pstruct_out.volume /
            pstruct_out.composition.get_reduced_composition_and_factor()[1],
            'nform':
            pstruct_out.composition.get_reduced_composition_and_factor()[1],
            'spg_in': spg_in,
            'spg_out': spg_out,
            'spg_in_high_tol': spg_in_high_tol,
            'spg_out_high_tol': spg_out_high_tol,
            **inpd
        })
    return pd.DataFrame(records)


def get_relax_records(groups: orm.Group,
                      encut: int,
                      group_with_structures=True,
                      **extra_filters) -> pd.DataFrame:
    """
    Query a group of structure and gather valid relaxation results

    Args:
        group (list): A list of group to query from
        encut (float): The cut off energy
        group_with_structures (bool): Wether the group contains structures or relaxation workchains directly.
    
    Any extra keyword arguments are applied as filters
    """

    q = QueryBuilder()
    q.append(
        Group,
        filters={'id': {
            'in': [g.id for g in groups]
        }},  # Groups from which the records are selected
        project=['label'],
    )
    if group_with_structures:
        q.append(
            Node, with_group=Group
        )  # The Node in the Group that is used as the input of the relaxation
        q.append(orm.WorkflowNode,
                 filters={'attributes.exit_status': {
                     'in': [0, 600]
                 }, 
                 'process_type': 'aiida.workflows:vaspu.relax'
                 },
                 project=['label', 'uuid'])  # Relaxaton workchain
    else:
        # The VaspRelaxation is directly in the group
        q.append(orm.WorkflowNode,
                 filters={'attributes.exit_status': {
                     'in': [0, 600]
                 },
                 'process_type': 'aiida.workflows:vaspu.relax'
                 },
                 with_group=Group,
                 project=['label', 'uuid'])  # Relaxaton workchain
    q.append(StructureData, with_outgoing=orm.WorkflowNode,
             project=['*', 'label'])  # The structure
    q.append(Str,
             with_outgoing=orm.WorkflowNode,
             edge_filters={'label': {
                 'like': '%potential_family'
             }},
             project=['attributes.value'])  # PP
    q.append(Dict,
             with_outgoing=orm.WorkflowNode,
             edge_filters={'label': {
                 'like': '%parameters'
             }},
             project=['attributes'],
             filters={
                 'or': [{
                     'attributes.vasp.encut': encut
                 }, {
                     'attributes.incar.encut': encut
                 }]
             })  # Support for both new and old style

    # Additional filters for the input parameters - support for both vasp and incar
    # overriding namespace
    for key, value in extra_filters.items():
        q.append(Dict,
                 with_outgoing=orm.WorkflowNode,
                 edge_filters={'label': {
                     'like': '%parameters'
                 }},
                 filters={
                     'or': [{
                         f'attributes.vasp.{key}': value
                     }, {
                         f'attributes.incar.{key}': value
                     }]
                 })

    q.append(Dict,
             with_incoming=orm.WorkflowNode,
             edge_filters={'label': 'misc'},
             project=['attributes.total_energies', 'attributes.maximum_force'
                      ])  # support for both new and old style
    q.append(StructureData, with_incoming=orm.WorkflowNode,
             edge_filters={'label':'relax__structure'},
             project=['*'])  # output structure for the relax
    ncount = q.count()
    print("Entries: {}".format(ncount))

    records = []
    for group_label, relax_label, relax_uuid, struct, label, pp, inpd, eng, fmax, relaxed_structure in tqdm(
            q.iterall(), total=ncount):
        # Deconstruct as some of the structures does not have vasp in the input
        if 'vasp' in inpd:
            inpd = inpd['vasp']
        else:
            inpd = inpd['incar']

        # Unpack the energy
        if 'energy_no_entropy' in eng:
            eng = eng['energy_no_entropy']
        else:
            eng = eng['energy_extrapolated']
        xc = get_functional(inpd, pp)

        pstruct_in = get_pmg_struct(struct)  # This is the input structure
        pstruct_out = get_pmg_struct(relaxed_structure)

        spg_in, _ = pstruct_in.get_space_group_info(1e-5, -1)
        spg_out, _ = pstruct_out.get_space_group_info(1e-5, -1)
        spg_in_high_tol, _ = pstruct_in.get_space_group_info(1e-1, -1)
        spg_out_high_tol, _ = pstruct_out.get_space_group_info(1e-1, -1)

        records.append({
            'group_label':
            group_label,
            'original_label':
            label,
            'original_uuid':
            struct.uuid,
            'structure':
            struct,
            'structure_uuid':
            struct.uuid,
            'label':
            label,
            'relax_label':
            relax_label,
            'relax_uuid':
            relax_uuid,
            'energy':
            eng,
            'functional':
            xc,
            'pp':
            pp,
            'max_force':
            fmax,
            'umap':
            get_u_map(struct, inpd.get('ldauu')),
            'note':
            'relax',
            'relax_structure':
            relaxed_structure,
            'pmg_struct':
            pstruct_in,
            'pmg_struct_relaxed':
            pstruct_out,
            'formula':
            pstruct_in.composition.reduced_formula,
            'volume':
            pstruct_out.volume,
            'volume':
            pstruct_out.volume,
            'volume_per_fu':
            pstruct_out.volume /
            pstruct_out.composition.get_reduced_composition_and_factor()[1],
            'nform':
            pstruct_out.composition.get_reduced_composition_and_factor()[1],
            'spg_in': spg_in,
            'spg_out': spg_out,
            'spg_in_high_tol': spg_in_high_tol,
            'spg_out_high_tol': spg_out_high_tol,
            **inpd
        })
    return pd.DataFrame(records)


def get_enumerate_records(groups: List[orm.Group], encut: int,
                          **extra_filters) -> pd.DataFrame:
    """
    Add the results from magnetic enumeration workchains
    Take the enumerations with the topography 
    Group -> StructureData -> Relax -> *SpinEnumeration*
    """
    # Query to included magnetic enumeration results
    q = QueryBuilder()
    q.append(Group,
             filters={'id': {
                 'in': [g.id for g in groups]
             }},
             project=['label'])
    q.append(Node, with_group=Group, project=['*', 'label'])
    # First relaxation
    q.append(orm.WorkflowNode, filters={
        'process_type': 'aiida.workflows:vaspu.relax'
    })
    q.append(StructureData)
    # Relaxed structure used for spin enumeration
    q.append(orm.WorkflowNode,
             filters={'attributes.exit_status': {
                 'in': [0]
             }, 'process_type': {'like': '%magnetic'}},
             tag='enum')
    q.append(orm.WorkflowNode,
             filters={'attributes.exit_status': {
                 'in': [0, 600]
             }, 
            'process_type': 'aiida.workflows:vaspu.relax'
             },
             project=['label', 'uuid'],
             tag='relax')
    q.append(StructureData, with_outgoing='relax',
             project=['*', 'label'])  # Input structure for the relax
    q.append(Str,
             with_outgoing='relax',
             edge_filters={'label': {
                 'like': '%potential_family'
             }},
             project=['attributes.value'])
    q.append(Dict,
             with_outgoing='relax',
             edge_filters={'label': {
                 'like': '%parameters'
             }},
             project=['attributes'],
             filters={
                 'or': [{
                     'attributes.vasp.encut': encut
                 }, {
                     'attributes.incar.encut': encut
                 }]
             })  # Support for both new and old style
    q.append(Dict,
             with_incoming=orm.WorkflowNode,
             edge_filters={'label': 'misc'},
             project=['attributes.total_energies', 'attributes.maximum_force'])
    q.append(StructureData, with_incoming='relax',
             project=['*'])  # output structure for the relax
    ncount = q.count()

    print("Entries for magnetic enumeration calculations: {}".format(ncount))

    records = []
    for group_label, struct, label, relax_label, relax_uuid, actual_struct, acutal_label, pp, inpd, eng, maximum_force, relaxed_structure in tqdm(
            q.iterall(), total=ncount):
        # Deconstruct as some of the structures does not have vasp in the input
        if 'vasp' in inpd:
            inpd = inpd['vasp']
        else:
            inpd = inpd['incar']

        # Unpack the energy
        if 'energy_no_entropy' in eng:
            eng = eng['energy_no_entropy']
        else:
            eng = eng['energy_extrapolated']
        xc = get_functional(inpd, pp)

        pstruct_in = get_pmg_struct(actual_struct)
        pstruct_out = get_pmg_struct(relaxed_structure)

        records.append({
            'group_label':
            group_label,
            'original_label':
            label,
            'original_uuid':
            struct.uuid,
            'structure':
            actual_struct,
            'structure_uuid':
            actual_struct.uuid,
            'label':
            label,
            'relax_label':
            relax_label,
            'relax_uuid':
            relax_uuid,
            'energy':
            eng,
            'pp':
            pp,
            'functional':
            'xc',
            'max_force':
            maximum_force,
            'umap':
            get_u_map(struct, inpd.get('ldauu')),
            'note':
            'relax',
            #'final_magnetization': tot_mag,
            'relax_structure':
            relaxed_structure,
            'pmg_struct':
            pstruct_in,
            'pmg_struct_relaxed':
            pstruct_out,
            'formula':
            pstruct_in.composition.reduced_formula,
            'volume':
            pstruct_out.volume,
            'volume_per_fu':
            pstruct_out.volume /
            pstruct_out.composition.get_reduced_composition_and_factor()[1],
            'nform':
            pstruct_out.composition.get_reduced_composition_and_factor()[1],
            **inpd
        })
    return pd.DataFrame(records)


# Standard list of the column and their definitions. The other columns are typically the
# INCAR parameters used in the calculation.
STD_COLUMN_DEFS = {
    'group_label': 'Label of the Group containing the structure',
    'original_label':
    'Original label of the structure, before any transformation',
    'original_uuid': 'The UUID of the original structure',
    'structure': 'StructureData instance of the calculation',
    'structure_uuid':
    'The UUID of the StructureData goes into the calculation',
    'label': 'Label of the final calculation that gives the energy',
    'relax_label': 'Label of the relaxation workchain',
    'relax_uuid': 'UUID of the relaxation workchain',
    'energy': 'Raw energy output',
    'functional': 'Name of the exchange-correlation functional',
    'pp': 'Pseudopotential family',
    'max_force': 'Maximum force in the relaxed structure',
    'umap': 'Mapping of the Hubbard U for each element',
    'note': 'Custom remark of the row',
    'relax_structure': 'Relaxed StructureData instance',
    'pmg_struct': 'Pymatgen Structure of the input StructureData',
    'formula': 'The reduced formula returned by pymatgen',
    'volume': 'Volume of the relaxed structure',
    'volume_per_fu': 'Volume per formula unit of the relaxed structure',
    'nform': 'Number of formula units'
}

MANDATORAY_COLUMNS = [
    'label', 'energy', 'formula', 'pmg_struct', 'pmg_struct_relaxed'
]
