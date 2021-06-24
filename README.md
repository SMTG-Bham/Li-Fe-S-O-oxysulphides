## Data for paper "Predicting lithium iron oxysulphides for battery cathodes"

Content of the repository

- `search-data` contains the relaxed structures sampled by ab initio random structure searching (AIRSS).
- `hull-structure-files` contains the structures further relaxed using VASP and used for the convex hull construction.
- `calcs` contains the raw calculation data export from AiiDA file repository for easy access. These data are also in the `all-calc.aiida` archive file.
- `all-calcs-data.aiida` is an AiiDA export file containing all calculation performed and reported in paper with their full provenance. This file can be imported into an AiiDA instance (`aiida-core >= 1.6.3`). Raw data may also be access by simply unzipping the file.
- `notebooks` contains the notebooks for data analysis and visualisation.
- `figures` contains the figures in the PNG and the SVG format. The SVG figures are generated using Inkscape.
- `visualisation` contains the structure files used for making the figure.


## Importing data into a new AiiDA Database

The provenance data stored in the `all-calc-data.aiida` can be imported into a AiiDA database and allow the analysis in the `notebooks` folder to be reproduced.
A few dependencies needs to be installed:

- `aiida-core>=1.6.3` is the mean package of the AiiDA framework.
- `aiida-vasp>2.0.1` contains the interface to VASP.
- `ase` and `pymatgen` are needed for the analysis. `ase==3.21` and `pymatgen==2020.10.20` were used in this project, but other compatible version should work as well.
- additional requires in the `requirements.txt`.

A few more steps are needed to step up a new AiiDA profile before importing the archive. Please follow the installation instructions on the [offical documentation](https://aiida.readthedocs.io/).

The magic "%aiida" can be enabled by a file `aiida_magic_register.py` in `<home_folder>/.ipython/profile_default/startup`:


```python
if __name__ == '__main__':

    try:
        import aiida
        del aiida
    except ImportError:
        # AiiDA is not installed in this Python environment
        pass
    else:
        from aiida.tools.ipython.ipython_magics import register_ipython_extension
        register_ipython_extension()
```

Note that you may want to run `ipython` once so the profile folers are created. If not, just 

Once everything is in place, the following commands will import all data:

```base
verdi archive import -G li-fe-s-o-paper -- all-calc-data.aiida
```

As each calculation/workflow is identified based on a unique ID (uuid), the analysis code in the notebook can be reused. 
