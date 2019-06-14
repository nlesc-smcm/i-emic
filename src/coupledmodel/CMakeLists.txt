set(COUPLEDMODEL_INCLUDE_DIRS
  ../ocean/
  ../atmosphere/
  ../seaice/
  ../utils/
  ../dependencygrid/
  ${CMAKE_CURRENT_SOURCE_DIR}
  )

add_library(coupledmodel SHARED CoupledModel.C)

target_compile_definitions(coupledmodel PUBLIC ${COMP_IDENT})

target_include_directories(coupledmodel PUBLIC ${COUPLEDMODEL_INCLUDE_DIRS})

target_link_libraries(coupledmodel PRIVATE
    ${MPI_CXX_LIBRARIES}
    ${Belos_LIBRARIES}
    ${Belos_TPL_LIBRARIES}
    ${Epetra_LIBRARIES}
    ${Epetra_TPL_LIBRARIES}
    ${ML_LIBRARIES}
    ${ML_TPL_LIBRARIES}
    atmosphere
    ocean
    seaice
    ifpack_mrilu
    trios
)

install(TARGETS coupledmodel DESTINATION lib)
