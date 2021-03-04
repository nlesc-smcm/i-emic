# Trying to find Trilinos
if (DEFINED ENV{TRILINOS_DIR} AND NOT DEFINED Trilinos_DIR)
  set (Trilinos_DIR $ENV{TRILINOS_DIR}/lib/cmake/Trilinos)
elseif (APPLE AND NOT DEFINED ENV{TRILINOS_DIR} AND NOT DEFINED Trilinos_DIR)
  # By default look in the homebrew install directory
  set (Trilinos_DIR /usr/local/opt/trilinos/lib/cmake/Trilinos/)
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
	"-g -O3 -Wall -ffixed-line-length-132 -fdefault-real-8 -fPIC -ffree-line-length-none")
  set (CMAKE_CXX_FLAGS "-g -O3 -Wall -fPIC -Wno-deprecated-declarations -DDEBUGGING_NEW")
  set (COMP_IDENT GNU)
  
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Intel") #ifort, icc
  
  set (CMAKE_Fortran_FLAGS "-g -r8 -O3 -fPIC -warn -DASCII_TOPO -DWITH_UNION -heap-arrays 1 -extend-source 132")
  set (CMAKE_CXX_FLAGS "-O3 -fPIC")
  set (COMP_IDENT INTEL)
  
endif ()

include_directories(SYSTEM ${Trilinos_INCLUDE_DIRS})
include_directories(SYSTEM ${Trilinos_TPL_INCLUDE_DIRS})

set(library_directories ${Trilinos_LIBRARY_DIRS})
list(APPEND library_directories ${Trilinos_TPL_LIBRARY_DIRS})

list(APPEND library_directories ${MRILU_DIR}/lib/) # MRILU

link_directories(${library_directories})

set(include_flags "")
foreach(dir IN LISTS Trilinos_INCLUDE_DIRS Trilinos_TPL_INCLUDE_DIRS)
    string(CONCAT include_flags "${include_flags}" " -I${dir}")
endforeach(dir)

set(link_flags "")
foreach(dir IN LISTS Trilinos_LIBRARY_DIRS Trilinos_TPL_LIBRARY_DIRS)
    string(CONCAT link_flags "${link_flags}" " -L${dir}")
endforeach()
foreach(lib IN LISTS Teuchos_LIBRARIES Teuchos_TPL_LIBRARIES Epetra_LIBRARIES
        Epetra_TPL_LIBRARIES EpetraExt_LIBRARIES EpetraExt_TPL_LIBRARIES MPI_CXX_LIBRARIES)
    if ("${lib}" MATCHES "^/")
        string(CONCAT link_flags "${link_flags}" " ${lib}")
    else ()
        string(CONCAT link_flags "${link_flags}" " -l${lib}")
    endif ()
endforeach()

file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/Config.mk "CXXFLAGS+=${include_flags}")
file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/Config.mk "\nLDFLAGS+=${link_flags}")
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/Config.mk DESTINATION include)
unset(include_flags)

file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/null.cpp "")
add_library(iemic SHARED ${CMAKE_CURRENT_BINARY_DIR}/null.cpp)
target_link_libraries(iemic PRIVATE
    ${MPI_CXX_LIBRARIES}
    ${Belos_LIBRARIES}
    ${Belos_TPL_LIBRARIES}
    ${Epetra_LIBRARIES}
    ${Epetra_TPL_LIBRARIES}
    ${ML_LIBRARIES}
    ${ML_TPL_LIBRARIES}
    trios
    ocean
)
set(IEMIC_INTERFACE
    coupledmodel
    globaldefs
    transient
    utils
)

find_package(JDQZPP)
if (JDQZPP_FOUND)
  include_directories(${JDQZPP_INCLUDE_DIRS})
  list(APPEND library_directories ${JDQZPP_LIBRARY_DIRS})
  list(APPEND IEMIC_INTERFACE ${JDQZPP_LIBRARIES})
  add_definitions(-DHAVE_JDQZPP)
endif ()

find_package(RAILS)

target_link_libraries(iemic PUBLIC continuation ${IEMIC_INTERFACE})
if (APPLE)
    list(TRANSFORM IEMIC_INTERFACE PREPEND "LINKER:-reexport_library,$<TARGET_FILE:")
    list(TRANSFORM IEMIC_INTERFACE APPEND ">")
    target_link_options(iemic PUBLIC ${IEMIC_INTERFACE})
endif ()
install(TARGETS iemic DESTINATION lib)

# ------------------------------------------------------------------
# Get targets and dependencies of the components
add_subdirectory(mrilucpp)

add_subdirectory(trios)
add_subdirectory(ocean)
add_subdirectory(atmosphere)
add_subdirectory(seaice)
add_subdirectory(coupledmodel)
add_subdirectory(transient)
add_subdirectory(continuation)

add_subdirectory(globaldefs)
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
