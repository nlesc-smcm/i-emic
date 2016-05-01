cmake_minimum_required(VERSION 2.8.12.2)
  
set(TESTS_INCLUDE_DIRS
  ../utils/
  ../ocean/
  ../atmosphere/
  ../coupledmodel/
  ../supervector/
  ../idrsolver/
  ../multiarray/
  ../gmressolver/
  ../thetastepper/
  ../newton/
  ../continuation/
  ${CMAKE_CURRENT_SOURCE_DIR}
  )

add_executable(test_supervec test_supervec.C)

target_compile_definitions(test_supervec PUBLIC ${COMP_IDENT})

target_include_directories(test_supervec PUBLIC ${TESTS_INCLUDE_DIRS})

target_link_libraries(test_supervec PUBLIC ${I-EMIC_LIBS})
target_link_libraries(test_supervec PUBLIC ${library_dependencies})

install(TARGETS test_supervec DESTINATION test/default)