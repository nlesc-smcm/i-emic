set(TEST_INCLUDE_DIRS
  ../topo/
  ../lyapunov/
  ../transient/
  ${CMAKE_CURRENT_SOURCE_DIR}
  )

set(TEST_SOURCES
  test_map.C
  test_domain.C
  test_vector.C
  test_topo.C
  test_ocean.C
  trns_ocean.C
  intt_ocean.C
  reft_ocean.C
  test_parameterlist.C
  intt_2dmoc.C
  test_atmos.C
  test_oceanatmos.C
  test_seaice.C
  test_coupled.C
  intt_coupled.C
  test_integrals.C
  test_matrix.C
  test_ams.C
  )

include(BuildExternalProject)
BuildExternalProject(
  GTest
  GIT_REPOSITORY git://github.com/google/googletest.git
  GIT_TAG release-1.8.0)

BuildExternalProject_find_package(GTest)

include_directories(BEFORE ${GTEST_INCLUDE_DIRS})

set(TEST_LIBRARIES
  iemic
  ${MPI_CXX_LIBRARIES}
  ${Belos_LIBRARIES}
  ${Belos_TPL_LIBRARIES}
  ${Epetra_LIBRARIES}
  ${Epetra_TPL_LIBRARIES}
  ${EpetraExt_LIBRARIES}
  ${EpetraExt_TPL_LIBRARIES}
  ${GTEST_BOTH_LIBRARIES}
  ${CMAKE_THREAD_LIBS_INIT}
)

if (RAILS_FOUND)
  list(APPEND TEST_LIBRARIES ${RAILS_LIBRARIES})
  list(APPEND TEST_INCLUDE_DIRS ${RAILS_INCLUDE_DIRS})
  add_definitions(-DHAVE_RAILS)
endif()

if (JDQZPP_found)
  list(APPEND TEST_SOURCES
    test_jdqz.C
    )
endif()

foreach(test_source ${TEST_SOURCES})
  get_filename_component(test_name ${test_source} NAME_WE)
  add_executable(${test_name} ${test_source} TestDefinitions.C)

  string(REPLACE "_" ";" test_name_list ${test_name})
  list(GET test_name_list 1 test_base_name)
  file(MAKE_DIRECTORY ${PROJECT_SOURCE_DIR}/test/${test_base_name})
  add_test(NAME ${test_name} COMMAND ${test_name} WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/test/${test_base_name})

  target_compile_definitions(${test_name} PUBLIC ${COMP_IDENT})
  target_include_directories(${test_name} PUBLIC ${TEST_INCLUDE_DIRS})
  target_link_libraries(${test_name} PUBLIC ${TEST_LIBRARIES})
endforeach()

execute_process(COMMAND ${MPIEXEC} --oversubscribe echo RESULT_VARIABLE status OUTPUT_QUIET)
if(NOT status)
    set(MPI_OVERSUBSCRIBE "--oversubscribe")
endif()

# additional parallel tests
get_filename_component(test_name test_matrix.C NAME_WE)
add_test(NAME partest_matrix_1 COMMAND ${MPIEXEC} -np 1 ${MPI_OVERSUBSCRIBE} ${CMAKE_CURRENT_BINARY_DIR}/${test_name}
  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/test/matrix)
add_test(NAME partest_matrix_2 COMMAND ${MPIEXEC} -np 2 ${MPI_OVERSUBSCRIBE} ${CMAKE_CURRENT_BINARY_DIR}/${test_name}
  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/test/matrix)
add_test(NAME partest_matrix_4 COMMAND ${MPIEXEC} -np 4 ${MPI_OVERSUBSCRIBE} ${CMAKE_CURRENT_BINARY_DIR}/${test_name}
  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/test/matrix)
add_test(NAME partest_matrix_8 COMMAND ${MPIEXEC} -np 8 ${MPI_OVERSUBSCRIBE} ${CMAKE_CURRENT_BINARY_DIR}/${test_name}
  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/test/matrix)

