cmake_minimum_required(VERSION 3.3)
  
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

target_link_libraries(coupledmodel PRIVATE ${MPI_CXX_LIBRARIES})

install(TARGETS coupledmodel DESTINATION lib)
