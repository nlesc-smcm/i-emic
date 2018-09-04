cmake_minimum_required(VERSION 2.8.12.2)

set(UTILS_INCLUDE_DIRS
  ${CMAKE_CURRENT_SOURCE_DIR}
  ../ocean/
  )

add_library(utils SHARED Utils.C)

target_link_libraries(utils PUBLIC ${library_dependencies})

target_include_directories(utils PUBLIC ${UTILS_INCLUDE_DIRS})

target_compile_definitions(utils PUBLIC ${COMP_IDENT})


install(TARGETS utils DESTINATION lib)
