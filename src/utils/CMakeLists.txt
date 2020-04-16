add_library(utils SHARED Utils.C Combined_MultiVec.C Model.C)

target_link_libraries(utils PRIVATE
    ${MPI_CXX_LIBRARIES}
    ${Teuchos_LIBRARIES}
    ${Teuchos_TPL_LIBRARIES}
    ${Epetra_LIBRARIES}
    ${Epetra_TPL_LIBRARIES}
    ${EpetraExt_LIBRARIES}
    ${EpetraExt_TPL_LIBRARIES}
)

target_link_libraries(utils PUBLIC globaldefs trios)

target_compile_definitions(utils PUBLIC ${COMP_IDENT})
target_include_directories(utils PUBLIC .)

install(FILES ComplexVector.H JDQZInterface.H Model.H Utils.H DESTINATION include)
install(TARGETS utils DESTINATION lib)
