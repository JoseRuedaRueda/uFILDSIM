# FILDSIM code

This code was created as a synthetic INPA, although it was enhanced to have FILDSIM capabilities. The code is written in fortran, although a complete set of python routines to prepare the inputs and process the output is written.

Up to now, these python libraries are distributed together with the ScintillatorSuite (<https://github.com/JoseRuedaRueda/ScintSuite>). This ScintSuite will be open sourced soon

## Citation of the code
To cite the code please use the FILDSIM [![DOI](https://zenodo.org/badge/758880137.svg)](https://zenodo.org/doi/10.5281/zenodo.12512086) plus the FILDSIM or INPASIM article, depending if you are running the code for FILD or for INPA:
- FILDSIM mode: J Galdon-Quiroga etal 2018 PlasmaPhys.Control.Fusion 60 105005
- INPASIM mode: J Rueda-Rueda et al 2024 Plasma Phys. Control. Fusion 66 065025

## Installation and documentation
### Cloning the code
Just go to the parent folder and clone the code. Note: When developing it, this upgraded FILDSIM code was the Synthetic INPA (SINPA) code. This is how it is called in the ScintSutie code. So better cloned in a SINPA folder, to keep the old convention of naming.

First, you need to say to git to do not record the file rights, as from the shares (IPP file system) to your personal laptop or the configuration in other machines, they can be different, so we can get fake file change right. Just open a terminal and write:
```bash
git config --global core.fileMode false
git config core.fileMode false
```

After it, you can manually and easily clone the code:

For AUG users: `<ParentFolderForSuite> = /shares/departments/AUG/users/$USER` (it can be changed if you want)
```bash
cd <ParentFolderForSuite>
git clone https://github.com/JoseRuedaRueda/uFILDSIM SINPA
cd SINPA
git checkout <branch>
```
### Compilation
In other to compile the code, you should use the version of gfortran included in gcc version 9.3.0 (GCC) [This is the default when you do module load gcc in toki, as of 4/1/2021]. Notice that any GCC between 8 and 10 seems to work.
There are 2 ways of compiling the code:
  - `make all`: It will compile the code in a way that is compatible with all processors, should work from your laptop to clusters
  - `make all_cluster`: It will compile it with some optimization flags though for high performance CPU and clusters. Notice that these options are cluster dependent, for example this compiled version work in MPCDF cluster toks but not toki. Special compilation for other clusters could be added upon request. (thos option does no longer work as it contained flags for the inter CPU, but toki was upgraded to AMD hardware)

### Test run
SINPA comes with a test input configuration file. Once your code is compiled, go to the file `./runs/InitialTest/inputs/InitialTest.cfg` and change the variables indicated with a !!! in the comments (they are the ones pointing to speficit paths of your installation, runDir and GeomDir). Once that is done, open the terminal and type:
``` bash
  cd <SINPA_FOLDER>
  ./bin/SINPA.go ./runs/InitialTest/inputs/InitialTest.cfg
```

The code should run. You will notice that there are some 'wrong markers' These are just markers which not collide with anything, for examples those with pitch too small (in degrees) are likely to don't collide with the scintillator. Don't worry, everything is fine

### Documentation
- All objects and methods are documented such that the user can understand what is going on, although the code is new and improvement in the documentation will come with new versions
- For examples of how to execute the code and analyze the outputs, see the examples python scripts in the ScintillatorSuite library (unfortunately, up to now this python library was not exported from the Suite and the complete ScintSute must be download and installed. In the future once the code development is frozen, they will be ported here)
- For a manual (in the process to be completed) with a complete benchmark with oldFILDSIM see <https://datashare.mpcdf.mpg.de/s/oQJw82fd705OEJ8>. Notice that this manual contain technical insformation about the models, inputs and outputs. It is not thought to be a guide to run the code and analyze the outputs

### Version control
Each release will be denoted by 2 numbers: a.b meaning:
    - b: bug fixed and improved comments and documentation. Some new capabilities could be added (see changelog). The higher the number, the better.
    - a: indicate major changes in the code, versions with different 'a' may be not compatible as input/output formats may have changed

## Develop
### Branches
- master: Stable branch, things should work, may be a delay including new features
- dev-branch: developers branch, may have some small bugs due to recent features being added. Include the latest features, not recommended for general public
- 'tmp'-branch: linked to specific commits to include new features or to keep track of not fully developed features. Do not use these branches except you are the developer in charge of the new feature. Unicorns can appear

### Note for developers
- Before changing anything in a module open a issue in GitLab to start a discussion

### Issues and new implementations
If a new implementation is required, open the appropriate issue in the GIT and link it to the milestone if it corresponds (if possible). The following tags are available:

- Documentation: improve the documentation of a given section.
- Feature request: request to implement a new feature in the code.
- Minor mod.: request to implement minor modifications in the code.
- Enhancement: modify the implementation of a given feature to improve the efficiency or make easier some processing.
- Discussion: a forum to discuss ideas of implementation.
- Bug: minor error found in the code. To be corrected at the earliest convenience.
- Major error: an important error has to be solved in the code as soon as possible.
- Minor priority: Label for maintainer, indicates that the request has low priority in the ToDo list

## Useful links
- FILDSIM code: <https://gitlab.mpcdf.mpg.de/jgq/FILDSIM.git>
- i-HIBPSIM code: <https://gitlab.mpcdf.mpg.de/poyo/ihibpsim>
- ScintillatorSuite library: <https://github.com/JoseRuedaRueda/ScintSuite>
