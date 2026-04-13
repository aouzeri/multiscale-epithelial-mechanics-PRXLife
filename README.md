# 3D Active gel vertex model for simulating inflating and deflating spherical epithelial domes


The present work is a [hiperlife](https://doi.org/10.5281/zenodo.14927572) [[1](#References)] based project for simulating inflating and deflating epithelial domes using vertex modelling. In this model [[2,3](#References)], the mechanical behavior of each cellular surface is described with an active gel model of the actomyosin cortex capturing its viscoelasticity, active contractility and turnover. Cell surfaces are triangulated resulting in a 3D vertex model with curved cellular surfaces and the model is solved on each of these surfaces using finite elements. To limit cell interpenetration, we included a repulsive potential between cells. Additionally, surfaces enclose fixed cellular volumes. 

## License 
Distributed under the GNU GENERAL PUBLIC LICENSE. See [LICENSE](LICENSE) for details.

## Hiperlife [[1](#References)]

We use a parallel finite element library called [hiperlife](https://doi.org/10.5281/zenodo.14927572) (High Performance Library for Finite Elements) [[1](#References)], which serves as the core numerical engine for the simulations presented in this study. The library is openly distributed to the community and aims at providing a computational framework to address problems of cell and tissue mechanobiology for a wide range of cases and users, with special focus on curved surfaces. The library is written in C++, uses the Message Passing Interface (MPI) paradigm for parallelism, and is built on top of several packages of the [Trilinos Project](https://trilinos.github.io/). The installation of the hiperlife libraries can be carried out by following the README.md and INSTALL.md files provided at https://doi.org/10.5281/zenodo.14927572. 

## Installation
See [INSTALL.md](INSTALL.md) for detailed installation instructions.

## Code Organization and Project Setup

The folder is organized as follows :
```bash
.                        # <--- top-level (or root) directory of the project (called agvm_domes)   
├── build                # <--- Temporary build system files (CMake, object files, etc.)
├── simulation           # <--- Folder containing the meshes to run a simulation
|   ├── results          # <--- Folder where to run the simulation
|   └── mesh             # <--- Folder containing the mesh files
├── src                  # <--- Folder containing the source code (implementation files)
|   ├── CMakeLists.txt   # <--- Application configuration file.
|   ├── AuxVXDomes.cpp   # <--- Implementation of helper functions
|   ├── AuxVXDomes.h     # <--- Header file declaring helper functions
|   └── VXDomes.cpp      # <--- Main, entry point of the program
├── post-processing      # <--- Contains the post-processing file used in Paraview (.pvsm)
├── CMakeLists.txt       # <--- Top-level configuration file.
├── userConfig.cmake     # <--- User's configuration file.
├── README.md            # <--- Project overview and usage instruction
├── INSTALL.md           # <--- Instructions on how to install and build the project 
├── AUTHORS.md           # <--- List of project authors and contributors  
└── .gitignore           # <--- Contains files to be ignored by git 
```

## Running a simulation

Simulations are run in the folder [simulation/results](simulation/results). To run a simulation, first copy all files contained in mesh/ to the folder results/

```bash
cd simulation/results
cp ../mesh/*.* .
```

The folder *simulation/mesh/* contains mesh files called vertexmesh_X.txt where X corresponds to the mesh number, as well as a configuration file [VMconfig.cfg](VMconfig.cfg) which stores all the necessary parameter values for the simulation.

To run a simulation :

```bash
mpirun -n 8 hlVXDomes VMconfig.cfg > output.out &
```

where *8* corresponds to the number of processors used, *hlVXDomes* the name of the executable, *> out8procs.out* to store the simulation output in a file.

IMPORTANT: if the path the executable is not explicitely specified, then it should be specified in ~/.bashrc, .e.g., by adding the following line to the file
```bash
export PATH="$HOME/.local/agvm_domes/:$PATH"
```
and the file needs to be sourced by
```bash
source ~/.bashrc
```
See [INSTALL.md](INSTALL.md) for more details.


## Data visualiwation and post-processing

At each timestep X, the code outputs a solutionX.vtm file and one solutionX.Y.vtu file per processor Y used which stores the values for the different degrees of freedom. For example, with 4 processors, step 89 outputs solution89.vtm and solution89.0.vtu, solution89.1.vtu, solution89.2.vtu, solution89.3.vtu, solution89.4.vtu. To visualize the results we use [Paraview](https://www.paraview.org/) which is an open source post-processing visualization engine freely available for download. 
We used version 5.11.2 available [here](https://www.paraview.org/download/) and we provide a [state file](post-processing/state.pvsm) which loads the solution*.vtm files and provides a cross section view of the the dome at different timesteps. 

To do so :

1) Open Paraview, e.g., using the executable located in *<path/to/paraview/bin/*
2) Load a state with *File > Load State*
3) Choose state.pvsm in the *post-processing* folder
4) In the dropdown menu, select *Choose File Names*
5) At the bottom of the Load State Options window (XMLMultiBlockDataReader), select the *solution..vtm* files to add the collection of solution files.

## Expected outputs

Here we show snapshots of the results for an inflated and deflating dome. We simulate spherical dome inflation and deflation in a monolayer attached to a substrate by fixing nodes outside a circular basal footprint. Once the dome deflates, nodes in renewed contact with the substrate are fixed.

We first start with a flat monolayer

![plot](post-processing/initial_state.png)

We then inflate and hold the dome over a large periode of time allowing sufficient time for the active gel dynamics to relax. 

![plot](post-processing/inflated_dome.png)

Finaly, we rapidely deflate the dome, causing the tissue to buckle and to form wrinkling patterns upon contact with the substrate. 

![plot](post-processing/deflating_dome_1.png)
![plot](post-processing/deflating_dome_2.png)
![plot](post-processing/deflating_dome_3.png)

Time shown in the snapshots are in turnover timescale. 

Note that the snapshots shown correspond to approximately 0h, 1h35,  2h, 3h15 and 10h of simulation time, respectively, using 8 cores on a Dell Intel® Xeon(R) Silver 4208 CPU @ 2.10GHz × 32. During the deflation phase, the time-step decreases for numerical convergence as the dome buckles and cells come in contact. 

## Cite 

You can acknowledge this package in any publication that relies on it using reference [[2]](#References).

## References

[1] Santos-Oliván, D., Vilanova, G., & Torres-Sánchez, A. (2025). hiperlife. Zenodo. https://doi.org/10.5281/zenodo.14927572

[2] Nimesh Chahare, Adam Ouzeri, Thomas Wilson, Pradeep K. Bal, Tom Golde, Guillermo Vilanova, Pau Pujol-Vives, Pere Roca-Cusachs, Xavier Trepat, Marino Arroyo. Multiscale wrinkling dynamics in epithelial shells. bioRxiv 2025.06.30.662426; doi: https://doi.org/10.1101/2025.06.30.662426

[3] Adam Ouzeri. Theory and computation of multiscale epithelial mechanics : from active gels to vertex models. Doctoral thesis. doi: https://doi.org/10.5821/dissertation-2117-402162

