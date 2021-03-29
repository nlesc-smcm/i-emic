file(COPY . DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/mrilu)

find_package(X11)
get_filename_component(X11_LIB_DIR "${X11_X11_LIB}" DIRECTORY)

message ("-- X11 include: " ${X11_INCLUDE_DIR})
message ("-- X11 library: " ${X11_LIB_DIR})

enable_language(Fortran)

if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(CMAKE_Fortran_FLAGS
    "-g -O3 -fPIC -ffree-line-length-none")
  set(MODFLAG "-J")
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(CMAKE_Fortran_FLAGS
    "-g -u -O3 -fPIC -DWITH_UNION -heap-arrays 1")
  set(MODFLAG "-module ")
endif ()

# XCode compilers require explicit sysroot to find C include headers
if (APPLE AND NOT CMAKE_C_COMPILER MATCHES "^/usr/bin/")
    set (SYSROOT_FLAG "-isysroot ${CMAKE_OSX_SYSROOT}")
endif ()

set(MRILU_DIR ${CMAKE_CURRENT_BINARY_DIR}/build PARENT_SCOPE)
# Turns out that cmake sometimes returns multiple dirs for X11
set(X11_INCLUDE_DIRS ${X11_INCLUDE_DIR})
list(TRANSFORM X11_INCLUDE_DIRS PREPEND "-I")
list(JOIN X11_INCLUDE_DIRS " " X11_INCLUDE_DIRS)

set(MAKEFILE_INC_CONTENT "
ARCH=
NUMBASE=${CMAKE_CURRENT_BINARY_DIR}/build
SOURCEDIR=${CMAKE_CURRENT_BINARY_DIR}
OSNAME=
MAKE=make
LD=${CMAKE_Fortran_COMPILER}
LDFLAGS=
CC=${CMAKE_C_COMPILER}
CFLAGS=-ansi ${SYSROOT_FLAG}
CPPFLAGS=${X11_INCLUDE_DIRS}
CXX=${CMAKE_CXX_COMPILER}
MPI_CC=
F77=${CMAKE_Fortran_COMPILER}
FPPFLAGS= 
FFLAGS=${CMAKE_Fortran_FLAGS}
MPI_F77=
F90=${CMAKE_Fortran_COMPILER}
MPI_F90=
MODFLAG=${MODFLAG}
LIB_X11=${X11_LIB_DIR}
ARFLAGS=-crsv")

if (NOT APPLE)
    set(MAKEFILE_INC_CONTENT "${MAKEFILE_INC_CONTENT}U")
endif ()

file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/mrilu/makefile.inc "${MAKEFILE_INC_CONTENT}")

add_custom_target(mrilu COMMAND make && make install
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mrilu)
