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
  ../topo/
  ${CMAKE_CURRENT_SOURCE_DIR}
  )

add_executable(test_jdqz     test_jdqz.C)
add_executable(test_topo     test_topo.C)
add_executable(test_ocean    test_ocean.C)
add_executable(test_coupled  test_coupled.C)
add_executable(test_atmos    test_atmos.C)


target_compile_definitions(test_jdqz PUBLIC ${COMP_IDENT})
target_compile_definitions(test_topo PUBLIC ${COMP_IDENT})
target_compile_definitions(test_ocean PUBLIC ${COMP_IDENT})
target_compile_definitions(test_coupled PUBLIC ${COMP_IDENT})
target_compile_definitions(test_atmos PUBLIC ${COMP_IDENT})


target_include_directories(test_jdqz PUBLIC ${TESTS_INCLUDE_DIRS})
target_include_directories(test_topo PUBLIC ${TESTS_INCLUDE_DIRS})
target_include_directories(test_ocean PUBLIC ${TESTS_INCLUDE_DIRS})
target_include_directories(test_coupled PUBLIC ${TESTS_INCLUDE_DIRS})
target_include_directories(test_atmos PUBLIC ${TESTS_INCLUDE_DIRS})

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


target_link_libraries(test_jdqz PUBLIC ${I-EMIC_LIBS})
target_link_libraries(test_jdqz PUBLIC ${library_dependencies})
target_link_libraries(test_jdqz PUBLIC ${GTEST_LIBRARIES})

target_link_libraries(test_topo PUBLIC ${I-EMIC_LIBS})
target_link_libraries(test_topo PUBLIC ${library_dependencies})
target_link_libraries(test_topo PUBLIC ${GTEST_LIBRARIES})

target_link_libraries(test_ocean PUBLIC ${I-EMIC_LIBS})
target_link_libraries(test_ocean PUBLIC ${library_dependencies})
target_link_libraries(test_ocean PUBLIC ${GTEST_LIBRARIES})

target_link_libraries(test_coupled PUBLIC ${I-EMIC_LIBS})
target_link_libraries(test_coupled PUBLIC ${library_dependencies})
target_link_libraries(test_coupled PUBLIC ${GTEST_LIBRARIES})

target_link_libraries(test_atmos PUBLIC ${I-EMIC_LIBS})
target_link_libraries(test_atmos PUBLIC ${library_dependencies})
target_link_libraries(test_atmos PUBLIC ${GTEST_LIBRARIES})

# Pthread
find_package(Threads)
target_link_libraries(test_jdqz     PUBLIC ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(test_topo     PUBLIC ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(test_ocean    PUBLIC ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(test_coupled  PUBLIC ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(test_atmos    PUBLIC ${CMAKE_THREAD_LIBS_INIT})

install(TARGETS test_jdqz     DESTINATION test/jdqz)
install(TARGETS test_topo     DESTINATION test/topo)
install(TARGETS test_ocean    DESTINATION test/ocean)
install(TARGETS test_coupled  DESTINATION test/coupled)
install(TARGETS test_atmos    DESTINATION test/atmos)
