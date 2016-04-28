cmake_minimum_required(VERSION 2.8.12.2)
project(i-emic)

enable_language (Fortran)

# Trilinos root directory should be in PATH
find_package(Trilinos)

set(CMAKE_CXX_COMPILER ${Trilinos_CXX_COMPILER})
set(CMAKE_Fortran_COMPILER ${Trilinos_Fortran_COMPILER})

include_directories(${Trilinos_INCLUDE_DIRS})
include_directories(${Trilinos_TPL_INCLUDE_DIRS})

set(library_dependencies ${Trilinos_LIBRARIES})
list(APPEND library_dependencies ${Trilinos_TPL_LIBRARIES})

set(library_directories ${Trilinos_LIBRARY_DIRS})
list(APPEND library_directories ${Trilinos_TPL_LIBRARY_DIRS})

add_subdirectory(ocean)