cmake_minimum_required(VERSION 2.8.12.2)
  
set(SUPERVECTOR_INCLUDE_DIRS
  ../utils/
  ../ocean/
  ${CMAKE_CURRENT_SOURCE_DIR}
  )

add_library(supervector SuperVector.C)

target_compile_definitions(supervector PUBLIC ${COMP_IDENT})

target_include_directories(supervector PUBLIC ${SUPERVECTOR_INCLUDE_DIRS})
