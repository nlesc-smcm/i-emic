cmake_minimum_required(VERSION 3.3)
project(i-emic)

set(DATA_DIR "'${CMAKE_CURRENT_SOURCE_DIR}/data/'")
message("-- Data is stored in: " ${DATA_DIR})

# MRILU_DIR should be available
if (NOT DEFINED MRILU_DIR)
  if (NOT DEFINED ENV{MRILU_DIR})
	message(FATAL_ERROR "Define MRILU_DIR environment variable such that we can find <MRILU_DIR>/mod/, <MRILU_DIR>/lib/ etc. ")
  else ()
	SET(MRILU_DIR $ENV{MRILU_DIR})
  endif ()
endif ()
message("-- Looking for MRILU in ${MRILU_DIR}")

enable_testing()

add_subdirectory(cmake)
message("-- CMAKE module path: ${CMAKE_MODULE_PATH}")

add_subdirectory(src)

# If no CMAKE_INSTALl_PREFIX is available we set it to the current source dir
if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set (CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_SOURCE_DIR}/" CACHE PATH "default install path" FORCE )
endif()

message("-- Install prefix = ${CMAKE_INSTALL_PREFIX}")
