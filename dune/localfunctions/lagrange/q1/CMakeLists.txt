set(q1dir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/lagrange/q1)

set(q1_HEADERS  q1localbasis.hh q1localcoefficients.hh q1localinterpolation.hh)

# include not needed for CMake
# include $(top_srcdir)/am/global-rules

install(FILES ${q1_HEADERS} DESTINATION ${q1dir})
