set(MAIN_INCLUDE_DIRS
  ../continuation/
  ../topo/
  ../lyapunov/
  ../transient/
  ${CMAKE_CURRENT_SOURCE_DIR}
  )

set(MAIN_SOURCES
  run_ocean.C
  run_coupled.C
  time_ocean.C
  time_coupled.C
  run_topo.C
  run_ams.C
  )

set(MAIN_LIBRARIES
  iemic
  ${MPI_CXX_LIBRARIES}
  ${Belos_LIBRARIES}
  ${Belos_TPL_LIBRARIES}
  ${Epetra_LIBRARIES}
  ${Epetra_TPL_LIBRARIES}
  ${EpetraExt_LIBRARIES}
  ${EpetraExt_TPL_LIBRARIES}
)

if (RAILS_FOUND)
  list(APPEND MAIN_SOURCES run_lyapunov.C)
  list(APPEND MAIN_LIBRARIES ${RAILS_LIBRARIES})
  list(APPEND MAIN_INCLUDE_DIRS ${RAILS_INCLUDE_DIRS})
endif()

foreach(main_source ${MAIN_SOURCES})
  get_filename_component(name ${main_source} NAME_WE)
  add_executable(${name} ${main_source})
  target_link_libraries(${name} PUBLIC ${MAIN_LIBRARIES})
  target_include_directories(${name} PUBLIC ${MAIN_INCLUDE_DIRS})
  install(TARGETS ${name} DESTINATION bin)
endforeach()
