set(UTILS_INCLUDE_DIRS
  ${CMAKE_CURRENT_SOURCE_DIR}
  ../ocean/
  )

add_library(utils SHARED Utils.C GlobalDefinitions.C)

target_link_libraries(utils PRIVATE
    ${MPI_CXX_LIBRARIES}
    ${Epetra_LIBRARIES}
    ${Epetra_TPL_LIBRARIES}
    ${EpetraExt_LIBRARIES}
    ${EpetraExt_TPL_LIBRARIES}
)

target_include_directories(utils PUBLIC ${UTILS_INCLUDE_DIRS})

target_compile_definitions(utils PUBLIC ${COMP_IDENT})

install(TARGETS utils DESTINATION lib)
