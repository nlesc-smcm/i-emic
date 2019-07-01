add_library(globaldefs SHARED GlobalDefinitions.C)

target_link_libraries(globaldefs PRIVATE
    ${MPI_CXX_LIBRARIES}
    ${Teuchos_LIBRARIES}
    ${Teuchos_TPL_LIBRARIES}
    ${Epetra_LIBRARIES}
    ${Epetra_TPL_LIBRARIES}
)

target_compile_definitions(globaldefs PUBLIC ${COMP_IDENT})
target_include_directories(globaldefs PUBLIC .)

install(FILES GlobalDefinitions.H DESTINATION include)
install(TARGETS globaldefs DESTINATION lib)
