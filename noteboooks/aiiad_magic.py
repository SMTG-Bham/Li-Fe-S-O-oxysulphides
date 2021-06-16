try:
    import aiida
except ImportError:
    pass
else:
#    from toolchest import export_atoms
    import numpy as np
    import matplotlib.pyplot as plt
    import IPython

    if aiida.__version__.startswith('1'):
        from aiida.tools.ipython.ipython_magics import load_ipython_extension

        # Get the current Ipython session
        ipython = IPython.get_ipython()

        # Register the line magic
        load_ipython_extension(ipython)

try:
    import pymatgen
    import ase
except ImportError:
    pass
else:
    from pymatgen.io.ase import AseAtomsAdaptor
    from pymatgen.core import Structure
    from ase.visualize import view as aview

    def view(atoms, *args, **kwargs):
        """
        Allow viewing pymatgen.core.Structure using ase
        """
        if isinstance(atoms, (list, tuple)):
            if isinstance(atoms[0], Structure):
                atoms = [AseAtomsAdaptor.get_atoms(s) for s in atoms]
        elif isinstance(atoms, Structure):
            atoms = AseAtomsAdaptor.get_atoms(atoms)
        return aview(atoms, *args, **kwargs)


def extract_bash_alias():
    """
    Extract aliass from the bashrc file
    """
    import re
    import os
    alias_re = re.compile(r'^alias +(\w+)="(.+)"$')
    alias = []
    with open(os.environ["HOME"] + "/.bashrc") as fh:
        for line in fh:
            m = alias_re.match(line)
            if m:
                alias.append((m.group(1), m.group(2)))
    return alias


# Add alias define in bashrc to the
ipython = get_ipython()
bash_alias = extract_bash_alias()

for alias, cmd in bash_alias:
    ipython.magic("alias {} {}".format(alias, cmd))
    #eval("{}='{}'".format(alias, cmd))

del ipython
del bash_alias
