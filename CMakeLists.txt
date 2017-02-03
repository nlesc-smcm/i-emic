cmake_minimum_required(VERSION 2.8.12.2)
project(i-emic)

# SHARED_DIR should be available
if (DEFINED ENV{SHARED_DIR})
  set(DATA_DIR "'$ENV{SHARED_DIR}/i-emic/data/'")
else()
  message (FATAL_ERROR "SHARED_DIR environment variable should be available (absolute directory)")
endif()

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

add_subdirectory(src)

# If no CMAKE_INSTALl_PREFIX is available we set it to the current source dir
if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set (CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_SOURCE_DIR}/" CACHE PATH "default install path" FORCE )
endif()

message("-- Install prefix = ${CMAKE_INSTALL_PREFIX}")
