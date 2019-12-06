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
