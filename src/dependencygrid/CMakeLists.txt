cmake_minimum_required(VERSION 3.3)

set(DEPENDENCYGRID_INCLUDE_DIRS
  ${CMAKE_CURRENT_SOURCE_DIR}
  ../utils/
  )

add_library(dependencygrid SHARED DependencyGrid.C)

target_compile_definitions(dependencygrid PUBLIC ${COMP_IDENT})

install(TARGETS dependencygrid DESTINATION lib)
