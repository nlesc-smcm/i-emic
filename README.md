# Build instructions

## Dependencies:

- cmake
- lapack
- blas 
- openmpi       (`libopenmpi-dev` in ubuntu repository)
- hdf5-openmpi  (`libhdf5-openmpi-dev` in ubuntu repository)
- metis         (`wget http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/metis-5.1.0.tar.gz`)
- parmetis      (`wget http://glaros.dtc.umn.edu/gkhome/fetch/sw/parmetis/parmetis-4.0.3.tar.gz`)
- mrilu         (available in this repository)
- jdqzpp        (external project, fetched and installed by cmake)

## Environment variables
Define the following envorinment variables:

- `PLAT`: platform type, only used in a few batch scripts
- `SHARED_DIR`: base dir of project: `{SHARED_DIR}/i-emic`
- `MRILU_DIR`: installation directory of mrilu: `{MRILU_DIR/{lib,mod,...}` 


## Installation:
  - Install metis
  - Install parmetis (requires openmpi)
  - Install MRILU
  - Patch Trilinos to enable `Ifpack_MRILU`, see `notes/patch`
  - Install Trilinos, see example `*.cmake` files in `notes/`  
  In `*.cmake` file: adjust `METIS_LIBRARY_DIRS`,`TPL_METIS_INCLUDE_DIRS`
  `ParMETIS_LIBRARY_DIRS` and `TPL_ParMETIS_INCLUDE_DIRS`
	  


