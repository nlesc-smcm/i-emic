set(ATMOSPHERE_INCLUDE_DIRS
  ../ocean/
  ../seaice/
  ../utils/
  ../dependencygrid/
  ${CMAKE_CURRENT_SOURCE_DIR}
  )

add_library(atmosphere STATIC AtmosLocal.C Atmosphere.C)

target_compile_definitions(atmosphere PUBLIC ${COMP_IDENT})

target_link_libraries(atmosphere dependencygrid)
target_link_libraries(atmosphere ocean)
target_link_libraries(atmosphere seaice)

target_include_directories(atmosphere PUBLIC ${ATMOSPHERE_INCLUDE_DIRS})

install(TARGETS atmosphere DESTINATION lib)
