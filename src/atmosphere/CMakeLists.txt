cmake_minimum_required(VERSION 2.8.12.2)
  
set(ATMOSPHERE_INCLUDE_DIRS
  ../ocean/
  ../utils/
  ../multiarray/
  ../supervector/
  ../gmressolver/
  ${CMAKE_CURRENT_SOURCE_DIR}
  )

add_library(atmosphere Atmosphere.C)

target_compile_definitions(atmosphere PUBLIC ${COMP_IDENT})

target_include_directories(atmosphere PUBLIC ${ATMOSPHERE_INCLUDE_DIRS})
