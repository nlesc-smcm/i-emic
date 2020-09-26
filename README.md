# I-EMIC: an Implicit Earth system Model of Intermediate Complexity.
The I-EMIC is a parallel Earth system model that allows the use of large-scale dynamical systems techniques. The term  'implicit' refers to this model's main design feature: all model equations and constraints are solved simultaneously. A Jacobian matrix and preconditioning techniques are available for efficient Newton iterations that are involved in implicit time stepping and continuation strategies.

The I-EMIC contains a number of fully implicit submodels coupled through a modular framework. At the center of the coupled model is the implicit primitive equation ocean model THCM [1]. THCM is extended with models of modest complexity that describe varying aspects of the climate, including atmospheric heat and moisture transport, the formation of sea ice and the adjustment of albedo over snow and ice.

For a description of the model equations, see [2].

[1] Dijkstra, H. A., Oksuzoglu, H., Wubs, F. W., and Botta, E. F. F. (2001). A fully implicit model of the three-dimensional thermohaline ocean circulation. Journal of Computational Physics, 173(2):685–715.

[2] Mulder, T. E. Design and bifurcation analysis of implicit Earth System Models. Dissertation (June, 2019),
Utrecht University, https://dspace.library.uu.nl/handle/1874/381139

# Build instructions

## Dependencies:

|                |                                                                                            |
| -------------- | ------------------------------------------------------------------------------------------ |
| Cmake          | (version 2.8.12.2 or higher)                                                               |
| lapack         | (`liblapack-dev` in ubuntu repository, `mkl` on intel)                                     |
| blas           | (`libblas-dev` in ubuntu repository, `mkl` on intel)                                       |
| openmpi        | (`libopenmpi-dev` in ubuntu repository)                                                    |
| hdf5-openmpi   | (`libhdf5-openmpi-dev` in ubuntu repository, `hdf5/impi` on intel)                         |
| metis          | (`wget http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/metis-5.1.0.tar.gz`)                |
| parmetis       | (`wget http://glaros.dtc.umn.edu/gkhome/fetch/sw/parmetis/parmetis-4.0.3.tar.gz`)          |
| mrilu          | (available in this repository)                                                             |
| Trilinos       | <https://trilinos.org/download/>  (this project is tested to work with 11.12/11.14/12.12)  |
| jdqzpp         | <https://github.com/erik808/jdqzpp>                                                        |
| gtest          | (external project, fetched and installed by cmake)                                         |

### Compilers
Depending on architecture: ifort, gfortran, mpicc, mpicpc, mpic++, etc...

## Installation:

### Manually
  * Install cmake, lapack, blas, openmpi and hdf5-openmpi as described above.

  * Install metis

  * Install parmetis (requires openmpi)

  * Install Trilinos with support for metis and parmetis:
    * Create build directory `{trilinos_base_dir}/build`
    * Create cmake script `build/{something}.cmake`, for examples see `notes/trilinos_cmake_examples`
      * Adjust `METIS_LIBRARY_DIRS`, `TPL_METIS_INCLUDE_DIRS`, `ParMETIS_LIBRARY_DIRS` and `TPL_ParMETIS_INCLUDE_DIRS`.

    * Make cmake script executable and run it, install Trilinos
      * Possible failures: no lapack, blas or hdf5 libs. `hdf5-openmpi` might install in `/usr/include/hdf5/openmpi`, so you could extend `CPATH` and `LD_LIBRARY_PATH` appropriately: e.g.: `export CPATH=$CPATH:/usr/include/hdf5/openmpi` `export LIBRARY_PATH=$LIBRARY_PATH:/usr/lib/x86_64-linux-gnu/hdf5/openmpi`

  * Install JDQZPP

  * Install I-EMIC
    * Create build directory
    * Create cmake script, see for examples `notes/i-emic_cmake_examples`
    * Run cmake script
    * make install -j<#procs>

### Using EasyBuild
  * Make sure no PATHs are set in your `.bashrc`.

  * Clone the EasyBuild repository located at <https://github.com/nlesc-smcm/easybuild>

  * Install the Trilinos module. For Cartesius, a job script is provided.

  * In your `.bashrc`, load the newly installed modules and set the appropriate PATHs and compilers, e.g.
  ```
  module load 2019 CMake Trilinos/12.14.1-intel-2019a-Python-3.7.2

  export PATH=${HOME}/local/bin:${PATH}
  export LIBRARY_PATH=${HOME}/local/lib:${HOME}/local/lib64:${LIBRARY_PATH}
  export LD_LIBRARY_PATH=${HOME}/local/lib64:${HOME}/local/lib:${LD_LIBRARY_PATH}:$EBROOTIMKL/mkl/lib/intel64

  export CC=mpiicc
  export CXX=mpiicpc
  export FC=mpiifort
  ```

  * Install JDQZPP

  * Compile I-EMIC, e.g.
  ```
  mkdir build
  cd build
  cmake ~/i-emic
  make -j<#procs>
  ```

# General remarks
- See the test code for examples ^^
