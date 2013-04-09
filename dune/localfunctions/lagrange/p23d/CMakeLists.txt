set(p23ddir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/lagrange/p23d)

set(p23d_HEADERS  p23dlocalbasis.hh p23dlocalcoefficients.hh p23dlocalinterpolation.hh)

# include not needed for CMake
# include $(top_srcdir)/am/global-rules

install(FILES ${p23d_HEADERS} DESTINATION ${p23ddir})
