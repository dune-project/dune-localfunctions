set(prismp1dir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/lagrange/prismp1)

set(prismp1_HEADERS  prismp1localbasis.hh prismp1localcoefficients.hh prismp1localinterpolation.hh)

# include not needed for CMake
# include $(top_srcdir)/am/global-rules

install(FILES ${prismp1_HEADERS} DESTINATION ${prismp1dir})
