cmake_minimum_required(VERSION 2.8.12.2)
project(i-emic)

# Trying to find Trilinos
if (DEFINED ENV{TRILINOS_DIR} AND NOT DEFINED Trilinos_DIR)
  set (Trilinos_DIR $ENV{TRILINOS_DIR}/lib/cmake/Trilinos)
endif ()
message ("-- Trilinos directory: " ${Trilinos_DIR})

# ------------------------------------------------------------------
# Trilinos
find_package(Trilinos)

if (NOT Trilinos_FOUND)
  message(FATAL_ERROR " Trilinos root directory should be in PATH or TRILINOS_DIR.")
endif ()

set(CMAKE_CXX_COMPILER ${Trilinos_CXX_COMPILER})
set(CMAKE_Fortran_COMPILER ${Trilinos_Fortran_COMPILER})

enable_language(Fortran)

message("-- CMAKE_Fortran_COMPILER_ID ${CMAKE_Fortran_COMPILER_ID}")

if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU") #gfortran, gcc
  
  set (CMAKE_Fortran_FLAGS
	"-g -u -O0 -ffixed-line-length-132 -fdefault-real-8 -ffree-line-length-none -DASCII_TOPO")
  set (CMAKE_CXX_FLAGS "-g -std=c++11 -O0 -Wall -Wno-deprecated-declarations -DDEBUGGING_NEW")
  set (COMP_IDENT GNU)
  
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Intel") #ifort, icc
  
  set (CMAKE_Fortran_FLAGS "-g -r8 -O3 -DASCII_TOPO -DWITH_UNION -heap-arrays 1 -extend-source 132")
  set (CMAKE_CXX_FLAGS "-std=c++11 -O3")
  set (COMP_IDENT INTEL)
  
endif ()

include_directories(${Trilinos_INCLUDE_DIRS})
include_directories(${Trilinos_TPL_INCLUDE_DIRS})

set(library_directories ${Trilinos_LIBRARY_DIRS})
list(APPEND library_directories ${Trilinos_TPL_LIBRARY_DIRS})

list(APPEND library_directories ${MRILU_DIR}/lib/) # MRILU

list(APPEND library_dependencies ${Trilinos_LIBRARIES})
list(APPEND library_dependencies ${Trilinos_TPL_LIBRARIES})

# Locate JDQZPP external project
#  maintained as a git submodule but that is not used here (yet).
include(ExternalProject)
set_directory_properties(PROPERTIES EP_PREFIX ${CMAKE_SOURCE_DIR}/src/external)

ExternalProject_Add(jdqzpp
  GIT_REPOSITORY https://github.com/erik808/jdqzpp
  GIT_TAG origin/master
  PREFIX ${CMAKE_SOURCE_DIR}/src/external/jdqzpp
  )

ExternalProject_Get_Property(jdqzpp source_dir)
set(JDQZPP_INCLUDE_DIR ${source_dir}/include)
set(JDQZPP_LIBRARY_DIR ${source_dir}/lib)
message("-- JDQZ++ source dir: " ${source_dir})

include_directories(${JDQZPP_INCLUDE_DIR})
list(APPEND library_directories ${JDQZPP_LIBRARY_DIR})

link_directories(${library_directories})

# ------------------------------------------------------------------
# Additional libraries
set(EXTRA_LIBS
  ifpack_mrilu
  mrilucpp
  precon
  mtstor
  misc
  iosrc
  gfortran
  hdf5
  )

list(APPEND library_dependencies ${EXTRA_LIBS})

# ------------------------------------------------------------------
# Internal I-EMIC shared libraries
set(I-EMIC_LIBS
  coupledmodel
  atmosphere
  ocean
  supervector
  utils
  jdqz_tools
  )

# ------------------------------------------------------------------
# Get targets and dependencies of the components
add_subdirectory(ocean)	  
add_subdirectory(mrilucpp)
add_subdirectory(atmosphere)
add_subdirectory(coupledmodel)
add_subdirectory(supervector)
add_subdirectory(utils)
add_subdirectory(main)
add_subdirectory(tests)

# ------------------------------------------------------------------
# CCACHE
find_program(CCACHE_FOUND ccache)
if(CCACHE_FOUND)
  set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
  set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ccache)
  message("-- CCACHE found and applied")
else ()
  message("-- CCACHE not found!")
endif(CCACHE_FOUND)
