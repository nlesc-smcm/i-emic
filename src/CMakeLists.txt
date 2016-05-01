cmake_minimum_required(VERSION 2.8.12.2)
project(i-emic)

enable_language (Fortran)

find_package(Trilinos)

if (NOT Trilinos_FOUND)
  message(FATAL_ERROR " Trilinos root directory should be in PATH ")
endif ()

set(CMAKE_CXX_COMPILER ${Trilinos_CXX_COMPILER})
set(CMAKE_Fortran_COMPILER ${Trilinos_Fortran_COMPILER})

if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU") #gfortran, gcc
  set (CMAKE_Fortran_FLAGS
	"-g -u -O0 -ffixed-line-length-132 -fdefault-real-8 -ffree-line-length-none -DASCII_TOPO")
  set (CMAKE_CXX_FLAGS "-g -std=c++11 -O0 -Wall")
  set (COMP_IDENT GNU)  
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Intel") #ifort, icc
  set (CMAKE_Fortran_FLAGS "-g -r8 -O3 -DASCII_TOPO -DWITH_UNION -heap-arrays 1 -extend-source 132")
  set (CMAKE_CXX_FLAGS "-std=c++11 -O3 -DMPICH_SKIP_MPICXX")
  set (COMP_IDENT INTEL)
endif (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")

find_program(CCACHE_FOUND ccache)
if(CCACHE_FOUND)
  set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
  set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ccache)
endif(CCACHE_FOUND)

include_directories(${Trilinos_INCLUDE_DIRS})
include_directories(${Trilinos_TPL_INCLUDE_DIRS})

set(library_directories ${Trilinos_LIBRARY_DIRS})
list(APPEND library_directories ${Trilinos_TPL_LIBRARY_DIRS})
list(APPEND library_directories $ENV{MRILU_DIR}/lib/) # MRILU
list(APPEND library_directories $ENV{LD_LIBRARY_PATH}) 

link_directories(${library_directories})

set(library_dependencies ${Trilinos_LIBRARIES})
list(APPEND library_dependencies ${Trilinos_TPL_LIBRARIES})

set(EXTRA_LIBS
  hdf5
  ifpack_mrilu
  mrilucpp
  precon
  mtstor
  misc
  iosrc
  gfortran
  )

list(APPEND library_dependencies ${EXTRA_LIBS})

set(I-EMIC_LIBS
  ocean
  atmosphere
  coupledmodel
  supervector
  ifpack_mrilu
  mrilucpp
  utils
  )

add_subdirectory(ocean)
add_subdirectory(mrilucpp)
add_subdirectory(atmosphere)
add_subdirectory(coupledmodel)
add_subdirectory(supervector)
add_subdirectory(utils)
add_subdirectory(main)
add_subdirectory(tests)
