cmake_minimum_required(VERSION 3.3)

execute_process(COMMAND mkdir -p ${CMAKE_BINARY_DIR}/external/cmake/ TIMEOUT 5)
execute_process(COMMAND wget -O ${CMAKE_BINARY_DIR}/external/cmake/BuildExternalProject.cmake https://raw.githubusercontent.com/Sbte/BuildExternalProject/master/BuildExternalProject.cmake TIMEOUT 5)

if (NOT EXISTS "${CMAKE_BINARY_DIR}/external/cmake/BuildExternalProject.cmake")
  # Download BuildExternalProject Module
  file(DOWNLOAD
    https://raw.githubusercontent.com/Sbte/BuildExternalProject/master/BuildExternalProject.cmake
    ${CMAKE_BINARY_DIR}/external/cmake/BuildExternalProject.cmake
    TIMEOUT 5)
endif ()

set(CMAKE_MODULE_PATH "${CMAKE_BINARY_DIR}/external/cmake" PARENT_SCOPE)
