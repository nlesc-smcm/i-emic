if (APPLE)
    cmake_minimum_required(VERSION 3.13)
else ()
    cmake_minimum_required(VERSION 3.5)
endif ()
project(i-emic)
set(CMAKE_CXX_STANDARD 14)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(DATA_DIR "'${CMAKE_CURRENT_SOURCE_DIR}/data/'")
message("-- Data is stored in: " ${DATA_DIR})

enable_testing()

add_subdirectory(cmake)
message("-- CMAKE module path: ${CMAKE_MODULE_PATH}")

if (NOT DEFINED MRILU_DIR)
  if (NOT DEFINED ENV{MRILU_DIR})
    add_subdirectory(mrilu)
  else ()
    SET(MRILU_DIR $ENV{MRILU_DIR})
  endif ()
endif ()

find_program(CCACHE_FOUND ccache)
if(CCACHE_FOUND)
  set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
  set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ccache)
endif(CCACHE_FOUND)

add_subdirectory(src)

# If no CMAKE_INSTALl_PREFIX is available we set it to the current source dir
if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set (CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_SOURCE_DIR}/" CACHE PATH "default install path" FORCE )
endif()

message("-- Install prefix = ${CMAKE_INSTALL_PREFIX}")

set(CTEST_OPTIONS --force-new-ctest-process --verbose CACHE STRING "ctest options")

# Adding custom test target "check" because it is the only way to show the full output of test errors on console
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} ${CTEST_OPTIONS})
