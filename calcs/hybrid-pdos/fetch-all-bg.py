from aiida_vasp.parsers.file_parsers.vasprun import VasprunParser
from pathlib import Path

base = Path('.')
calc = base.glob('Li*')

for folder in calc:
    subcalcs = folder.glob('relax*')
    for subcalc in subcalcs:
        vasprun_path = subcalc / 'vasprun.xml'
        parser = VasprunParser(file_path=str(vasprun_path))
        print(subcalc, ':{:.2f} eV'.format(parser.band_properties['band_gap']))