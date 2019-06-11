set(SOURCES
  Transient.cpp
  ScoreFunctions.C)

add_library(ams SHARED ${SOURCES})

set(AMS_INCLUDE_DIRS
  ../utils/
  ${CMAKE_CURRENT_SOURCE_DIR}
  )

target_include_directories(ams PUBLIC ${AMS_INCLUDE_DIRS})

install(TARGETS ams DESTINATION lib)
