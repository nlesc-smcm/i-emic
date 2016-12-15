# Build instructions

## Dependencies:

- cmake
- lapack        (`liblapack-dev` in ubuntu repository, `mkl` on intel)
- blas          (`libblas-dev` in ubuntu repository, `mkl` on intel)
- openmpi       (`libopenmpi-dev` in ubuntu repository)
- hdf5-openmpi  (`libhdf5-openmpi-dev` in ubuntu repository, `hdf5/impi` on intel)
- metis         (`wget http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/metis-5.1.0.tar.gz`)
- parmetis      (`wget http://glaros.dtc.umn.edu/gkhome/fetch/sw/parmetis/parmetis-4.0.3.tar.gz`)
- mrilu         (available in this repository)
- Trilinos      <https://trilinos.org/download/>
- jdqzpp        (external project, fetched and installed by cmake)
- gtest

### Compilers
Depending on architecture: ifort, gfortran, mpicc, mpicpc, mpic++, etc... 

## Environment variables
Define the following environment variables:

- `PLAT`: platform type, only used in a few batch scripts
- `SHARED_DIR`: base dir of project: `{SHARED_DIR}/i-emic`
- `MRILU_DIR`: installation directory of mrilu: `{MRILU_DIR}/{lib,mod,...}` 


## Installation:
  * Install gtest

  * Install metis

  * Install parmetis (requires openmpi)

  * Install MRILU
	*  Create an include file `mrilu/makefile.inc`, for examples see `mrilu/makefile_inc_examples`
  
  * Patch Trilinos to enable `Ifpack_MRILU`, see `notes/trilinos_ifpackmrilu_patch`
	* For Trilinos versions 11.12 and 11.14 you can just copy `{version}/Ifpack.*` to `packages/ifpack/src` 	
	
  * Install Trilinos with support for metis and parmetis:
	* Create build directory `{trilinos_base_dir}/build`
	* Create cmake script `build/{something}.cmake`, for examples see `notes/trilinos_cmake_examples`

		* Adjust `METIS_LIBRARY_DIRS`, `TPL_METIS_INCLUDE_DIRS`, `ParMETIS_LIBRARY_DIRS` and `TPL_ParMETIS_INCLUDE_DIRS`.
  
		* Pass an extra flag to ifpack: `-D Ifpack_CXX_FLAGS:STRING="-DHAVE_IFPACK_MRILU" \`
		
		* Let `src/mrilucpp/Ifpack_MRILU.h` be available in the `PATH` when compiling Trilinos.

	* Make cmake script executable and run it, install Trilinos 
	  * Possible failures: no lapack, blas or hdf5 libs. `hdf5-openmpi` might install in `/usr/include/hdf5/openmpi`, so you could extend `CPATH` and `LD_LIBRARY_PATH` appropriately.


  * Install I-EMIC
	* Create build directory
	* Create cmake script
