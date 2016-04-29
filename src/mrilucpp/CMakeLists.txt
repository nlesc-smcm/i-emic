cmake_minimum_required(VERSION 2.8.12.2)
  
set(MRILUCPP_INCLUDE_DIRS
  ../utils/
  $ENV{MRILU_DIR}/mod
  ${CMAKE_CURRENT_SOURCE_DIR}
  )

add_library(mrilucpp mrilucpp.F90)
add_library(ifpack_mrilu Ifpack_MRILU.cpp)

target_compile_definitions(mrilucpp PUBLIC HAVE_IFPACK_MRILU)
target_compile_definitions(ifpack_mrilu PUBLIC HAVE_IFPACK_MRILU)

target_include_directories(ifpack_mrilu PUBLIC ${MRILUCPP_INCLUDE_DIRS})
target_include_directories(mrilucpp PUBLIC ${MRILUCPP_INCLUDE_DIRS})
