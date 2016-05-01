cmake_minimum_required(VERSION 2.8.12.2)
  
set(TESTS_INCLUDE_DIRS
  ../utils/
  ../ocean/
  ../atmosphere/
  ../coupledmodel/
  ${CMAKE_CURRENT_SOURCE_DIR}
  )

add_executable(test_supervec test_supervec.C)

target_compile_definitions(test_supervec PUBLIC ${COMP_IDENT})

target_include_directories(test_supervec PUBLIC ${TESTS_INCLUDE_DIRS})
