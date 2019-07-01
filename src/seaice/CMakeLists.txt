add_library(seaice STATIC SeaIce.C)

target_compile_definitions(seaice PUBLIC ${COMP_IDENT})

target_include_directories(seaice PUBLIC .)

target_link_libraries(seaice PRIVATE atmosphere dependencygrid ocean utils)
