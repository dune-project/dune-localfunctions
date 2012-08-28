set(raviartthomas2q2ddir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/raviartthomas/raviartthomas2q2d)

set(raviartthomas2q2d_HEADERS      raviartthomas2q2dlocalbasis.hh
                                raviartthomas2q2dlocalcoefficients.hh
                                raviartthomas2q2dlocalinterpolation.hh)

# include not needed for CMake
# include $(top_srcdir)/am/global-rules

install(FILES ${raviartthomas2q2d_HEADERS} DESTINATION ${raviartthomas2q2ddir})
