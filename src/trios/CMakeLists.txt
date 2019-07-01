add_library(trios STATIC
  TRIOS_Domain.C
  TRIOS_BlockPreconditioner.C
  TRIOS_Saddlepoint.C
  TRIOS_SolverFactory.C
  TRIOS_Static.C
)

target_link_libraries(trios PRIVATE globaldefs ifpack_mrilu)
target_include_directories(trios PUBLIC .)
target_compile_definitions(trios PUBLIC ${COMP_IDENT})

install(FILES TRIOS_Domain.H DESTINATION include)
