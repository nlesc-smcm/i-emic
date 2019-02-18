# Build instructions

## Dependencies:

|                |                                                                                    |
| -------------- | ---------------------------------------------------------------------------------- |
| Cmake          | (version 2.8.12.2 or higher)                                                       |
| lapack         | (`liblapack-dev` in ubuntu repository, `mkl` on intel)                             |
| blas           | (`libblas-dev` in ubuntu repository, `mkl` on intel)                               |
| openmpi        | (`libopenmpi-dev` in ubuntu repository)                                            |
| hdf5-openmpi   | (`libhdf5-openmpi-dev` in ubuntu repository, `hdf5/impi` on intel)                 |
| metis          | (`wget http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/metis-5.1.0.tar.gz`)        |
| parmetis       | (`wget http://glaros.dtc.umn.edu/gkhome/fetch/sw/parmetis/parmetis-4.0.3.tar.gz`)  |
| mrilu          | (available in this repository)                                                     | 
| Trilinos       | <https://trilinos.org/download/>  (this project is tested to work with v11.12/14)  |
| jdqzpp         | (external project, fetched and installed by cmake)                                 |
| gtest          |  <https://github.com/google/googletest/>                                           |

### Compilers
Depending on architecture: ifort, gfortran, mpicc, mpicpc, mpic++, etc... 

## Installation:
  * Install cmake, lapack, blas, openmpi and hdf5-openmpi as described above.

  * Install gtest (ubuntu instructions, assuming sudo rights)
    *  `sudo apt-get install libgtest-dev`
    *  `cd /usr/src/gtest`
    *  `sudo cmake CMakeLists.txt`
    *  `sudo make`
    *  `sudo cp *.a /usr/lib`

  * Install metis

  * Install parmetis (requires openmpi)

  * Install MRILU
      *  Create an include file `mrilu/makefile.inc`, for examples see `mrilu/makefile_inc_examples`
      *  May depend on `libx11-dev` and `x11proto-core-dev`

  * Install Trilinos with support for metis and parmetis:
	* Create build directory `{trilinos_base_dir}/build`
	* Create cmake script `build/{something}.cmake`, for examples see `notes/trilinos_cmake_examples`

		* Adjust `METIS_LIBRARY_DIRS`, `TPL_METIS_INCLUDE_DIRS`, `ParMETIS_LIBRARY_DIRS` and `TPL_ParMETIS_INCLUDE_DIRS`.

	* Make cmake script executable and run it, install Trilinos 
	  * Possible failures: no lapack, blas or hdf5 libs. `hdf5-openmpi` might install in `/usr/include/hdf5/openmpi`, so you could extend `CPATH` and `LD_LIBRARY_PATH` appropriately: e.g.: `export CPATH=$CPATH:/usr/include/hdf5/openmpi` `export LIBRARY_PATH=$LIBRARY_PATH:/usr/lib/x86_64-linux-gnu/hdf5/openmpi`

  * Install I-EMIC
	* Create build directory
	* Create cmake script, see for examples `notes/i-emic_cmake_examples`
	* Run cmake script
    * make install -j<#procs>

# General remarks
- See the test code for examples ^^
