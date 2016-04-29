cmake_minimum_required(VERSION 2.8.12.2)

set(UTILS_INCLUDE_DIRS
  ${CMAKE_CURRENT_SOURCE_DIR}
  )

add_library(utils Utils.C)

target_compile_definitions(utils PUBLIC ${COMP_IDENT})
