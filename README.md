# Info
This is an example version of Molecular Dynamics of Expansion for GPU (MDEXPGPU for short). This is part of ongoing research of nucleation (the process of droplet formation) relevan in multiple areas such as weather forecast, steam turbine efficiency or covid pandemic. The software is custom built for execution on GPU for purpose of performing molecular simulation of SPC/E water model.
## Target problem
The intended application is to directly simulate the expansion. So in contrary to the traditionally emplyed technique of static system emplying methods controlling temperature(thermostats), which influence the quality of produced results, this approach aims to simulate the expanding volume closely via modelling the real situation like it happens in the supersonic expasnion apparature. We therefore simulate a parcel of expanding water vapor and observe the process of nucleation occuring within the simulated system.
## Demands
This presents a specific set of demands our software needs to satisfy.
 - custom iterative schema with expansion protocol adjsuting the volume of the system based on calculated position in the expanding nozzle
 - strict energy conservation
 	- this is required because no termostating is emplyed
 - long simulation times of <= 1 micro second
 	- this corresponds to 1 000 000 000 time steps (1 femto second time step)
 - larger number of particles (tested expansions for 1024 - 16384)
 	 - depending on the target droplet size and available hardware
 - fast execution (lower previous ~3 months of calculation)
 	- up to 108x speed up vs similar CPU implementation (optimized for CPU)

Demand for speed motivated the choice for GPGPU solution shown in this code showcase.

This is the older version of the code. More detailed research information are available from dissertation: *Development of Parallel Algorithms for Molecular Dynamics Simulation of Heterogeneous Atomistic Systems.* from David Celny, 2024. This work is available from library of FNSPE CTU [and can be downloaded here](hhttps://dspace.cvut.cz/handle/10467/114021).

This is part of my research project, if you want to extend its capabilities or find some issues with it, you can contact me [here](mailto:celny.david@gmail.com).

# Capabilities
 - NVE, NVT, EXPANSION ensambles/regimes
 	- set in during compilation
 - inlude BERENDSEN thermostating for testing purposes
 - include GEAR, MTK intagration schemes intended primarily for testing purposes
 	- not of much use unles you understand what you are doing
 - Argon, Nitrogen, SPC/E water model
 	- set in during compilation
 	- each has predefined interaction cutoff, verlet cutoff and elsectrostatic cutoff radius
 		 - Argon (3.0, 3.0\*1.5) Angstroms
 		 - Nitrogen (9.9234, 9.9234\*1.25) Angstroms
 		 - SPC/E (12.6622315608, 12.6875\*1.5, 12.6875) Angstroms
 - single or double precision controlled via macro switch
 	- higher precision has better energy conservation but is slower
 	- double precision is default (its better to wait bit more for more reliable results)
 - variable simualtion parameter settings
 	- Number of molecules (rounded up to nearest multiple of 16)
 	- initial density
 	- box size
 	- initial temperature estimate (temperature is consequence from velocities of molecules/atoms)
 	- integration parameters (like shake number, timestep time)
 	- holding structure settings (verlet list allocation size = maximum verlet size during simulation)
 	- total simulation time
 - multiple available output files generation
 	- energy, history, forces, velocity
 	- for more help on the file format please see [MACSIMUS](https://github.com/kolafaj/MACSIMUS/tree/main) manual
 		- plb, cp, cfg, cfa, log
 - initial configuration generation
 - supports loading of the configurations
 	- from cfg,cfa files (see [MACSIMUS](https://github.com/kolafaj/MACSIMUS/tree/main) manual)
 - supoorts for cuda ARCH 600

# How to instruction
If you want to compile and run this software on your machine you require Nvidia CUDA and at least CUDA 6.0 capable device.
The installation further assumes Linux operating system (This OS was used for development and testing - you may have some luck with making it run under windows but be aware paths convention is unix).

The following assumes that you have CUDA installed and enabled on your system. You can verify this by executing `nvidia-smi` (make sure you get no error and you can se your device in the list).
 Additionaly you require nvcc in compatible version with CUDA and at least version 17 ( you can verify this with 'nvcc --version' and CUDA documentation for compatibility).

## Compile
The project is compiled using CMake. This means you need the followng procedure.
 - Navigate root directory with the project (e.g. root CMakeList.txt is located)
 - create build directory if none is present `mkdir build`
 - configure cmake with `cmake -S. -B build`
 	- for specific configuration use the `cmake -B build/ -S. [-D<OPTION_NAME>=<OPTION_VALUE>]`
 		 - enabling float float with `cmake -B build/ -S. -DUSEFLOAT=ON`
		 - for different regime i.e. EXPANSION use `cmake -B build/ -S. -DSIM_TYPE_NAME=EXPANSION`
			- possible values inlude: CONSTANT, EXPANSION, BERENDSEN, MTK, GEAR
	- Dont worry about misspelling the SIM_TYPE_NAME you will get error with correct options if you make a mistake
		 - alternatively you can use the cmake-gui to configure - Be sure to delete from past changes when it does not work
 - compile the program for desired substance (default are all three substances:`argon`, `nitrogen`, `spce`)
 	 - from within the root of project call `cmake --build bild/ <name_of_substance>`
	 - for parallel compilation (recomended) use `-j<N_threads>` flag depending on your machine

## Execute
To execute you have to provide the program with file containing your parameters specification in form of `*.cdef` file.
Some of such files are provided in `examples` folder. The general structure of the files is in form of key=value pairs.
More specific example with expalanation is provided below:
```
N[0]=<number of molecules rounded up to 16>
rho=<density in kg/m3>
T=<temperature in K>
verlet_neighbour=<number of initial number of neighbors one molecule can have within the cutoff >
verlet_refresh=<how many timesteps should be performed before verlet list is refreshed, higher number means more risk some interactions are missed, usually 20 is good number>
dt=<timestep time in pico seconds usually 0.001 (enough for SPCE may be too fine for argon>
total_time=<total simulated time in pico seconds (start with 1 and see how it goes)>
n_shake=<number of iteration of shake bond conservation, 15 is good for SPCE, nitrogen does not need more than 2 and argon does not need any but 1 is fine to be specified)>
dt.eng=<specify the intervals for printout/update of this ouput file type in timestep (0 means only initial, -1 means last itration, and any positive number gives period of repetition)>

```
Program is then executed from build folder as `./gpu_module_CONSTANT_spce_double <path to your *.cdef>`
 - be aware that files of the same name as cdef will be created with the outputs requested in the folder with the cdef file
 - be also aware that based on the density of printout (lower frequency = higher density) quite large number of files may be generated or alternatively quite large files in size may be generated
 	- it safer to test more reasonable values first before you flood your disk with GB of text files
 - be mindfull that printing out information requires comminication and time which degrades the performance so printing every step is not good for the speed

 ## Execution printout
 - sucessful execution will print basic info about your gpu and loaded properties
 - if you see a waning with cuttoff than it means your system size is small copared to the verlet cutoff which may have physical detriment for the results generated - consult molecular dynamics books for more detailed information
 - the print is followed with debug information about verlet list size and mean interaction (mean number of neighbors) per molecule and potential realigns of this verlet
 	- this can give you more idea for choice for verlet_neighbour
 - when the simulation finishes you get the timing profile with total time and time spent on different operations
 	- there is also time per iteration which may help you for estimating how long (real time) it would take to simulate longer period (simulation time)
 	- you can test various verlet settings and sizes to see the influence

here is example how such result prinout can look like:
```
*** INFO - parsing parameters from : ../examples/SPCE/liquid_1024.cdef
*** 1024 molecules of SPCE water***
*** potential: 20 ,in FLOAT ***
*** BLOCK: 192, THREAD: 16, WORKLOAD: 64 ***
*** 31.289360*31.289360*31.289360 ***
*** cutoffs: norm=  12.662231, verlet=  19.031250, elstat=  12.687500 ***
*** V_width= 0, V_reset= 10, n_shake= 15 ***
*** t= 1.00e+00, dt= 1.0000e-03 ***
*** T = 300.00 ***
*** Rho = 1000.00 ***
*** 00:18:33 - 07.05.2024 ***
*** output: E=1, H=0, T=0, T2=0 ***
*** Legion5 ***
*** NVIDIA GeForce RTX 3070 Ti Laptop GPU ***
*** SM: 46, Gmem: 8192 MB, Smem/B: 48 kB, Cmem: 64 kB ***
!!! WARNING : The cutoff_verlet(19.031250) is higher then half of a box sizes(x/2: 15.644680, y/2: 15.644680, z/2: 15.644680). !!!
DEBUG ver_size: 32661, ver_size_new: 835888, ver_size_new_align: 835904, Int/mol: 816.312500
DEBUG ver_size: 0, ver_size_new: 835888, ver_size_new_align: 835904, Int/mol: 816.312500
DEBUG ver_size: 0, ver_size_new: 836276, ver_size_new_align: 836288, Int/mol: 816.687500
DEBUG ver_size: 0, ver_size_new: 837766, ver_size_new_align: 837824, Int/mol: 818.187500
DEBUG ver_size: 0, ver_size_new: 839134, ver_size_new_align: 839168, Int/mol: 819.500000
DEBUG ver_size: 0, ver_size_new: 840474, ver_size_new_align: 840512, Int/mol: 820.812500
DEBUG ver_size: 0, ver_size_new: 841752, ver_size_new_align: 841792, Int/mol: 822.062500
DEBUG ver_size: 0, ver_size_new: 842978, ver_size_new_align: 843008, Int/mol: 823.250000
DEBUG ver_size: 0, ver_size_new: 844126, ver_size_new_align: 844160, Int/mol: 824.375000
DEBUG ver_size: 0, ver_size_new: 845434, ver_size_new_align: 845440, Int/mol: 825.625000
=== TIMING PROFILE ===
t_total = 4.964727 s
t_setup = 2.116946 s, 42.639725 %
t_init  = 0.000236 s, 00.004754 %
t_iter  = 2.847386 s, 57.352318 %
 1_iter ~ 2.847386 ms
t_clean = 0.000159 s, 00.003203 %
```


# Concluding remarks
 This project was used throughout the work on dissertation and was originally designed for singular purpose of research of nucleation. During the history of the project from 2017 a lot has changed, which makes some portions of the code look and feel dated. Please excuse this issue of long runnig projects with very limited developer team.