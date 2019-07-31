add_library(dependencygrid STATIC DependencyGrid.C)

target_compile_definitions(dependencygrid PUBLIC ${COMP_IDENT})
target_include_directories(dependencygrid PUBLIC .)
