set(raviartthomas1q2ddir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/raviartthomas/raviartthomas1q2d)

set(raviartthomas1q2d_HEADERS      raviartthomas1q2dlocalbasis.hh
                                raviartthomas1q2dlocalinterpolation.hh
                                raviartthomas1q2dlocalcoefficients.hh)

# include not needed for CMake
# include $(top_srcdir)/am/global-rules

install(FILES ${raviartthomas1q2d_HEADERS} DESTINATION ${raviartthomas1q2ddir})
