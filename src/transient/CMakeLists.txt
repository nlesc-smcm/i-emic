set(SOURCES
  Transient.cpp
  ScoreFunctions.C)

add_library(transient SHARED ${SOURCES})

set(TRANSIENT_INCLUDE_DIRS
  ../utils/
  ${CMAKE_CURRENT_SOURCE_DIR}
  )

target_include_directories(transient PUBLIC ${TRANSIENT_INCLUDE_DIRS})

install(TARGETS transient DESTINATION lib)
