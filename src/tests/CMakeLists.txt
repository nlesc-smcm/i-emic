cmake_minimum_required(VERSION 2.8.12.2)
  
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
  test_domain.C
  test_vector.C
  test_jdqz.C
  test_topo.C
  test_ocean.C
  intt_ocean.C
  test_atmos.C
  test_seaice.C
  test_coupled.C
  intt_coupled.C
  )

# Google test
if (NOT GTEST_FOUND)
find_package(GTest)
endif ()

if (NOT GTEST_FOUND)
  if (DEFINED ENV{GTEST_DIR})
	message("-- GTest not automatically found, trying $ENV{GTEST_DIR}")
	set(GTEST_ROOT "$ENV{GTEST_DIR}")
  else ()
	message("-- GTest not automatically found, GTEST_DIR not set, trying $ENV{SHARED_DIR}")
	set(GTEST_ROOT "$ENV{SHARED_DIR}")
  endif ()
  
  unset(GTEST_INCLUDE_DIR CACHE) #FindGTEST quirk
  find_package(GTest)
endif ()

if (NOT GTEST_FOUND)
  message(FATAL_ERROR "GTest not found, tests will not compile")
endif ()

include_directories(${GTEST_INCLUDE_DIRS})



foreach(test_source ${TEST_SOURCES})
  get_filename_component(test_name ${test_source} NAME_WE)
  add_executable(${test_name} ${test_source})

  string(REPLACE "_" ";" test_name_list ${test_name})
  list(GET test_name_list 1 test_base_name)
  add_test(NAME ${test_name} COMMAND ${test_name} WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/test/${test_base_name})

  target_compile_definitions(${test_name} PUBLIC ${COMP_IDENT})
  target_include_directories(${test_name} PUBLIC ${TESTS_INCLUDE_DIRS})

  target_link_libraries(${test_name} PUBLIC ${I-EMIC_LIBS})
  target_link_libraries(${test_name} PUBLIC ${library_dependencies})
  target_link_libraries(${test_name} PUBLIC ${GTEST_BOTH_LIBRARIES})

  # Pthread
  target_link_libraries(${test_name} PUBLIC ${CMAKE_THREAD_LIBS_INIT})

endforeach()
