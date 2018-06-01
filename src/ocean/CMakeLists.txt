cmake_minimum_required(VERSION 2.8.12.2)

set(FORTRAN_SOURCES assemble.F90 atm.F90 ice.F90 monthly.F90
  boundary.F90 forcing.F90 grid.F90 inout.F90 itplbv.F90
  levitus.F90 mat.F90 matetc.F90 lev.F90 mix.F90
  res.F90 usr.F90 par.F90  global.F90 thcm_utils.F90
  scaling.F90 mix_imp.f mix_sup.F90 spf.F90 topo.F90
  usrc.F90 inserts.F90 probe.F90 integrals.F90)

set(CPP_SOURCES Ocean.C THCM.C OceanGrid.C
  TRIOS_Domain.C TRIOS_BlockPreconditioner.C TRIOS_Saddlepoint.C
  TRIOS_SolverFactory.C TRIOS_Static.C)

add_library(ocean STATIC ${FORTRAN_SOURCES} ${CPP_SOURCES})

set(OCEAN_INCLUDE_DIRS
  ../mrilucpp/
  ../utils/
  ../atmosphere/
  ../seaice/
  ../dependencygrid/
  ${CMAKE_CURRENT_SOURCE_DIR}
  )

target_include_directories(ocean PUBLIC ${OCEAN_INCLUDE_DIRS})
target_include_directories(ocean PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(ocean atmosphere)
target_link_libraries(ocean seaice)

target_compile_definitions(ocean PUBLIC DATA_DIR=${DATA_DIR} ${COMP_IDENT})

install(TARGETS ocean DESTINATION lib)
