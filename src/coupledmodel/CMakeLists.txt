add_library(coupledmodel SHARED CoupledModel.C)

target_compile_definitions(coupledmodel PUBLIC ${COMP_IDENT})

target_include_directories(coupledmodel PUBLIC .)

target_link_libraries(coupledmodel PUBLIC
    atmosphere
    ocean
    trios
    seaice
)

target_link_libraries(coupledmodel PRIVATE
    ${MPI_CXX_LIBRARIES}
    ${Belos_LIBRARIES}
    ${Belos_TPL_LIBRARIES}
    ${ML_LIBRARIES}
    ${ML_TPL_LIBRARIES}
    utils
)

install(TARGETS coupledmodel DESTINATION lib)
