add_library(atmosphere STATIC AtmosLocal.C Atmosphere.C)

target_compile_definitions(atmosphere PUBLIC ${COMP_IDENT})

target_link_libraries(atmosphere PRIVATE dependencygrid ocean seaice utils)

target_include_directories(atmosphere PUBLIC .)
