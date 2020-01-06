set(HEADERS
  bdfmcube.hh
)

add_subdirectory(cube)

install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/brezzidouglasfortinmarini)
