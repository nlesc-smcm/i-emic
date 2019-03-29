cmake_minimum_required(VERSION 3.3)
  
set(TESTS_INCLUDE_DIRS
  ../utils/
  ../ocean/
  ../atmosphere/
  ../seaice/
  ../coupledmodel/
  ../dependencygrid/
  ../continuation/
  ../topo/
  ${CMAKE_CURRENT_SOURCE_DIR}
  )

set(TEST_SOURCES
  test_map.C
  test_domain.C
  test_vector.C
  test_jdqz.C
  test_topo.C
  test_ocean.C
  intt_ocean.C
  reft_ocean.C
  test_atmos.C
  test_oceanatmos.C
  test_seaice.C
  test_coupled.C
  intt_coupled.C
  test_integrals.C
  test_matrix.C
  )

include(BuildExternalProject)
BuildExternalProject(
  GTest
  GIT_REPOSITORY git://github.com/google/googletest.git
  GIT_TAG release-1.8.0)

BuildExternalProject_find_package(GTest)

include_directories(BEFORE ${GTEST_INCLUDE_DIRS})

foreach(test_source ${TEST_SOURCES})
  get_filename_component(test_name ${test_source} NAME_WE)
  add_executable(${test_name} ${test_source})

  string(REPLACE "_" ";" test_name_list ${test_name})
  list(GET test_name_list 1 test_base_name)
  file(MAKE_DIRECTORY ${PROJECT_SOURCE_DIR}/test/${test_base_name})
  add_test(NAME ${test_name} COMMAND ${test_name} WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/test/${test_base_name})

  target_compile_definitions(${test_name} PUBLIC ${COMP_IDENT})
  target_include_directories(${test_name} PUBLIC ${TESTS_INCLUDE_DIRS})

  target_link_libraries(${test_name} PUBLIC ${I-EMIC_LIBS})
  target_link_libraries(${test_name} PUBLIC ${library_dependencies})
  target_link_libraries(${test_name} PUBLIC ${MPI_CXX_LIBRARIES})
  target_link_libraries(${test_name} PUBLIC ${GTEST_BOTH_LIBRARIES})

  # Pthread
  target_link_libraries(${test_name} PUBLIC ${CMAKE_THREAD_LIBS_INIT})

endforeach()

# additional parallel tests
get_filename_component(test_name test_matrix.C NAME_WE)
add_test(NAME partest_matrix_1 COMMAND mpirun -np 1 ${PROJECT_SOURCE_DIR}/build/src/tests/${test_name}
  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/test/matrix)
add_test(NAME partest_matrix_2 COMMAND mpirun -np 2 ${PROJECT_SOURCE_DIR}/build/src/tests/${test_name}
  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/test/matrix)
add_test(NAME partest_matrix_4 COMMAND mpirun -np 4 ${PROJECT_SOURCE_DIR}/build/src/tests/${test_name}
  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/test/matrix)
add_test(NAME partest_matrix_8 COMMAND mpirun -np 8 ${PROJECT_SOURCE_DIR}/build/src/tests/${test_name}
  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/test/matrix)

