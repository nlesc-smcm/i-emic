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
add_executable(test_jdqz test_jdqz.C)

target_compile_definitions(test_supervec PUBLIC ${COMP_IDENT})
target_compile_definitions(test_jdqz PUBLIC ${COMP_IDENT})

target_include_directories(test_supervec PUBLIC ${TESTS_INCLUDE_DIRS})
target_include_directories(test_jdqz PUBLIC ${TESTS_INCLUDE_DIRS})

# Google test
find_package(GTest)

if (NOT GTEST_FOUND)
  message("-- GTest not automatically found, trying $ENV{SHARED_DIR}")
  set(GTEST_ROOT "$ENV{SHARED_DIR}")
  unset(GTEST_INCLUDE_DIR CACHE) #FindGTEST quirk
  find_package(GTest)
endif ()
if (NOT GTEST_FOUND)
  message(FATAL_ERROR "GTest not found, tests will not compile")
endif ()

target_link_libraries(test_supervec PUBLIC ${I-EMIC_LIBS})
target_link_libraries(test_supervec PUBLIC ${library_dependencies})
target_link_libraries(test_supervec PUBLIC ${GTEST_LIBRARIES})

target_link_libraries(test_jdqz PUBLIC ${I-EMIC_LIBS})
target_link_libraries(test_jdqz PUBLIC ${library_dependencies})
target_link_libraries(test_jdqz PUBLIC ${GTEST_LIBRARIES})

# Pthread
find_package(Threads)
target_link_libraries(test_supervec PUBLIC ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(test_jdqz PUBLIC ${CMAKE_THREAD_LIBS_INIT})

install(TARGETS test_supervec DESTINATION test/default)
install(TARGETS test_jdqz DESTINATION test/jdqz)