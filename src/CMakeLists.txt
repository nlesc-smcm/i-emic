# Trying to find Trilinos
if (DEFINED ENV{TRILINOS_DIR} AND NOT DEFINED Trilinos_DIR)
  set (Trilinos_DIR $ENV{TRILINOS_DIR}/lib/cmake/Trilinos)
endif ()

# ------------------------------------------------------------------
# Trilinos
find_package(Trilinos REQUIRED)
message ("-- Trilinos directory: " ${Trilinos_DIR})

find_package(MPI REQUIRED)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})

enable_language(Fortran)

message("-- CMAKE_CXX_COMPILER:        ${CMAKE_CXX_COMPILER}")
message("-- CMAKE_CXX_COMPILER_ID:     ${CMAKE_CXX_COMPILER_ID}")
message("-- CMAKE_Fortran_COMPILER:    ${CMAKE_Fortran_COMPILER}")

if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU") #gfortran, gcc
  
  set (CMAKE_Fortran_FLAGS
	"-g -O0 -Wall -ffixed-line-length-132 -fdefault-real-8 -fPIC -ffree-line-length-none")
  set (CMAKE_CXX_FLAGS "-g -std=c++1y -O0 -Wall -fPIC -Wno-deprecated-declarations -DDEBUGGING_NEW")
  set (COMP_IDENT GNU)
  
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Intel") #ifort, icc
  
  set (CMAKE_Fortran_FLAGS "-g -r8 -O3 -fPIC -warn -DASCII_TOPO -DWITH_UNION -heap-arrays 1 -extend-source 132")
  set (CMAKE_CXX_FLAGS "-std=c++14 -O3 -fPIC")
  set (COMP_IDENT INTEL)
  
endif ()

include_directories(${Trilinos_INCLUDE_DIRS})
include_directories(${Trilinos_TPL_INCLUDE_DIRS})

set(library_directories ${Trilinos_LIBRARY_DIRS})
list(APPEND library_directories ${Trilinos_TPL_LIBRARY_DIRS})

list(APPEND library_directories ${MRILU_DIR}/lib/) # MRILU

find_package(JDQZPP REQUIRED)
if (JDQZPP_FOUND)
  include_directories(${JDQZPP_INCLUDE_DIR})
  list(APPEND library_directories ${JDQZPP_LIBRARY_DIR})
endif ()

link_directories(${library_directories})

# ------------------------------------------------------------------
# Internal I-EMIC libraries
set(I-EMIC_LIBS
  coupledmodel
  atmosphere
  seaice
  ocean
  utils
  dependencygrid
  jdqz_tools
  )

# ------------------------------------------------------------------
# Get targets and dependencies of the components
add_subdirectory(mrilucpp)

add_subdirectory(ocean)
add_subdirectory(atmosphere)
add_subdirectory(seaice)
add_subdirectory(coupledmodel)

add_subdirectory(utils)
add_subdirectory(dependencygrid)

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
