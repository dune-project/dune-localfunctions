set(HEADERS
  localbasis.hh
  localcoefficients.hh
  localinterpolation.hh
)

install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/brezzidouglasfortinmarini/cube)
