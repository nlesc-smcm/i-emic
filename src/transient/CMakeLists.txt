set(SOURCES
  Transient.cpp
  ScoreFunctions.C)

add_library(transient STATIC ${SOURCES})
target_include_directories(transient PUBLIC .)
target_link_libraries(transient PRIVATE
    ${MPI_CXX_LIBRARIES}
    ${Epetra_LIBRARIES}
    ${Epetra_TPL_LIBRARIES}
    ${EpetraExt_LIBRARIES}
    ${EpetraExt_TPL_LIBRARIES}
    ${Teuchos_LIBRARIES}
    ${Teuchos_TPL_LIBRARIES}
    utils
)
