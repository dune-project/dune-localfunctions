set(p0dir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/lagrange/p0)

set(p0_HEADERS  p0localbasis.hh p0localcoefficients.hh p0localinterpolation.hh)

# include not needed for CMake
# include $(top_srcdir)/am/global-rules

install(FILES ${p0_HEADERS} DESTINATION ${p0dir})
