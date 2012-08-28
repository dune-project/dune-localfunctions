set(pk2ddir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/lagrange/pk2d)

set(pk2d_HEADERS  pk2dlocalbasis.hh pk2dlocalcoefficients.hh pk2dlocalinterpolation.hh)

# include not needed for CMake
# include $(top_srcdir)/am/global-rules

install(FILES ${pk2d_HEADERS} DESTINATION ${pk2ddir})
