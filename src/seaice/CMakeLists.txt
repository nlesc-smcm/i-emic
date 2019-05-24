set(SEAICE_INCLUDE_DIRS
  ../ocean/
  ../atmosphere/
  ../utils/
  ../dependencygrid/
  ${CMAKE_CURRENT_SOURCE_DIR}
  )

add_library(seaice STATIC SeaIce.C)

target_compile_definitions(seaice PUBLIC ${COMP_IDENT})

target_link_libraries(seaice ocean)
target_link_libraries(seaice atmosphere)
target_link_libraries(seaice dependencygrid)

target_include_directories(seaice PUBLIC ${SEAICE_INCLUDE_DIRS})

install(TARGETS seaice DESTINATION lib)
