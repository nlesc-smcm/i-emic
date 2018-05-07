cmake_minimum_required(VERSION 2.8.12.2)
  
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

install(TARGETS coupledmodel DESTINATION lib)
