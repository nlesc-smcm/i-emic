cmake_minimum_required(VERSION 2.8.12.2)
project(i-emic)

if (DEFINED ENV{SHARED_DIR})
  set(DATA_DIR "'$ENV{SHARED_DIR}/i-emic/data/'")
else()
  message (FATAL_ERRORn "SHARED_DIR environment variable should be available")
endif()

message("-- Data is stored in: " ${DATA_DIR})

if (NOT DEFINED ENV{MRILU_DIR})
  message(FATAL_ERROR "Define MRILU_DIR environment variable such that we can find <MRILU_DIR>/mod/, <MRILU_DIR>/lib/ etc. ")
endif ()

add_subdirectory(src)
add_subdirectory(tests)

message("-- Executable will be put in " ${CMAKE_BINARY_DIR})