set(q2dir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/lagrange/q2)

set(q2_HEADERS  q2localbasis.hh q2localcoefficients.hh q2localinterpolation.hh)

# include not needed for CMake
# include $(top_srcdir)/am/global-rules

install(FILES ${q2_HEADERS} DESTINATION ${q2dir})
