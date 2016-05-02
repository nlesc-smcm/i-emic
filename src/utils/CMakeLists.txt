cmake_minimum_required(VERSION 2.8.12.2)

set(UTILS_INCLUDE_DIRS
  ${CMAKE_CURRENT_SOURCE_DIR}
  )

add_library(utils Utils.C)

set(FORTRAN_SOURCES
  ../jdqzpp/src/jdqz_tools/select.f
  ../jdqzpp/src/jdqz_tools/myexc.f
  )

add_library(jdqz_tools ${FORTRAN_SOURCES})
target_include_directories (jdqz_tools PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_compile_definitions(utils PUBLIC ${COMP_IDENT})
