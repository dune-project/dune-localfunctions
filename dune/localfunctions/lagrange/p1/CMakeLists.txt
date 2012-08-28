set(p1dir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/lagrange/p1)

set(p1_HEADERS  p1localbasis.hh p1localcoefficients.hh p1localinterpolation.hh)

# include not needed for CMake
# include $(top_srcdir)/am/global-rules

install(FILES ${p1_HEADERS} DESTINATION ${p1dir})
