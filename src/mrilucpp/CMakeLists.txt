set(MRILUCPP_INCLUDE_DIRS
  ../utils/
  ${MRILU_DIR}/mod
  ${CMAKE_CURRENT_SOURCE_DIR}
  )

add_library(mrilucpp SHARED mrilucpp.F90)
add_library(ifpack_mrilu SHARED Ifpack_MRILU.cpp)

target_compile_definitions(mrilucpp PUBLIC HAVE_IFPACK_MRILU)
target_compile_definitions(ifpack_mrilu PUBLIC HAVE_IFPACK_MRILU)

target_include_directories(ifpack_mrilu PUBLIC ${MRILUCPP_INCLUDE_DIRS})
target_include_directories(mrilucpp PUBLIC ${MRILUCPP_INCLUDE_DIRS})

link_directories(${MRILU_DIR}/lib)
list(APPEND FORTRAN_LIBS precon mtstor misc iosrc)
target_link_libraries(mrilucpp PRIVATE ${FORTRAN_LIBS})

target_link_libraries(ifpack_mrilu PRIVATE
    ${MPI_CXX_LIBRARIES}
    ${Epetra_LIBRARIES}
    ${Epetra_TPL_LIBRARIES}
    mrilucpp
    globaldefs
)
if (TARGET mrilu)
  add_dependencies(mrilucpp mrilu)
endif ()

install(TARGETS mrilucpp ifpack_mrilu DESTINATION lib)
