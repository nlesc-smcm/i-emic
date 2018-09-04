cmake_minimum_required(VERSION 3.3)

set(MAIN_INCLUDE_DIRS
  ../continuation/
  ../topo/
  ../ocean/
  ../atmosphere/
  ../seaice/
  ../mrilucpp/
  ../utils/
  ../dependencygrid/
  ${CMAKE_CURRENT_SOURCE_DIR}
  )

set(MAIN_SOURCES
  run_coupled.C
  run_ocean.C
  time_ocean.C
  run_topo.C
  )

foreach(main_source ${MAIN_SOURCES})
  get_filename_component(name ${main_source} NAME_WE)
  add_executable(${name} ${main_source})
  target_link_libraries(${name} PUBLIC ${I-EMIC_LIBS})
  target_link_libraries(${name} PUBLIC ${library_dependencies})
  target_include_directories(${name} PUBLIC ${MAIN_INCLUDE_DIRS})
  install(TARGETS ${name} DESTINATION bin)
endforeach()
