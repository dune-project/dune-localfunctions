set(qk_HEADERS
  qklocalbasis.hh
  qklocalcoefficients.hh
  qklocalinterpolation.hh)

install(FILES ${qk_HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/lagrange/qk)
