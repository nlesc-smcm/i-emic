set(FORTRAN_SOURCES assemble.F90 atm.F90 ice.F90
  boundary.F90 forcing.F90 grid.F90 inout.F90 itplbv.F90
  levitus.F90 mat.F90 matetc.F90 lev.F90 mix.F90
  res.F90 usr.F90 par.F90  global.F90 thcm_utils.F90
  scaling.F90 mix_imp.f mix_sup.F90 spf.F90 topo.F90
  usrc.F90 inserts.F90 probe.F90 integrals.F90)

set(CPP_SOURCES Ocean.C THCM.C OceanGrid.C)

add_library(ocean STATIC ${FORTRAN_SOURCES} ${CPP_SOURCES})
target_include_directories(ocean PUBLIC .)

target_link_libraries(ocean PRIVATE atmosphere ifpack_mrilu seaice utils)

target_compile_definitions(ocean PUBLIC DATA_DIR=${DATA_DIR} ${COMP_IDENT})

install(FILES Ocean.H DESTINATION include)
