cmake_minimum_required(VERSION 2.8.12.2)

add_executable(main main.C)

set(MAIN_INCLUDE_DIRS
  ${CMAKE_CURRENT_SOURCE_DIR}
  ../thetastepper/
  ../newton/
  ../continuation/
  )

target_link_libraries(main PUBLIC ${I-EMIC_LIBS})
target_link_libraries(main PUBLIC ${library_dependencies})
target_include_directories(main PUBLIC ${MAIN_INCLUDE_DIRS})
set_target_properties(main PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")