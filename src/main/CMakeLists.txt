cmake_minimum_required(VERSION 2.8.12.2)

add_executable(run_coupled run_coupled.C)
add_executable(run_ocean run_ocean.C)
add_executable(run_topo run_topo.C)


set(MAIN_INCLUDE_DIRS
  ${CMAKE_CURRENT_SOURCE_DIR}
  ../thetastepper/
  ../newton/
  ../continuation/
  ../topo/
  ../mrilucpp/
  )

target_link_libraries(run_coupled PUBLIC ${I-EMIC_LIBS})
target_link_libraries(run_coupled PUBLIC ${library_dependencies})
target_include_directories(run_coupled PUBLIC ${MAIN_INCLUDE_DIRS})

target_link_libraries(run_ocean PUBLIC ${I-EMIC_LIBS})
target_link_libraries(run_ocean PUBLIC ${library_dependencies})
target_include_directories(run_ocean PUBLIC ${MAIN_INCLUDE_DIRS})

target_link_libraries(run_topo PUBLIC ${I-EMIC_LIBS})
target_link_libraries(run_topo PUBLIC ${library_dependencies})
target_include_directories(run_topo PUBLIC ${MAIN_INCLUDE_DIRS})

install(TARGETS run_coupled DESTINATION bin)
install(TARGETS run_ocean DESTINATION bin)
install(TARGETS run_topo DESTINATION bin)
