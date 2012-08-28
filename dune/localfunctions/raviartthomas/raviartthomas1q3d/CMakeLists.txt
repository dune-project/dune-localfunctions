set(raviartthomas1q3ddir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/raviartthomas/raviartthomas1q3d)

set(raviartthomas1q3d_HEADERS      raviartthomas1q3dlocalbasis.hh
                                raviartthomas1q3dlocalcoefficients.hh
                                raviartthomas1q3dlocalinterpolation.hh)

# include not needed for CMake
# include $(top_srcdir)/am/global-rules

install(FILES ${raviartthomas1q3d_HEADERS} DESTINATION ${raviartthomas1q3ddir})
