set(pk1ddir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/lagrange/pk1d)

set(pk1d_HEADERS  pk1dlocalbasis.hh pk1dlocalcoefficients.hh pk1dlocalinterpolation.hh)

# include not needed for CMake
# include $(top_srcdir)/am/global-rules

install(FILES ${pk1d_HEADERS} DESTINATION ${pk1ddir})
