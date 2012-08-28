set(pk3ddir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/lagrange/pk3d)

set(pk3d_HEADERS  pk3dlocalbasis.hh pk3dlocalcoefficients.hh pk3dlocalinterpolation.hh)

# include not needed for CMake
# include $(top_srcdir)/am/global-rules

install(FILES ${pk3d_HEADERS} DESTINATION ${pk3ddir})
