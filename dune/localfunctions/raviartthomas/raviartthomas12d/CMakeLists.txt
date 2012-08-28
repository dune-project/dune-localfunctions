set(raviartthomas12ddir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/raviartthomas/raviartthomas12d)

set(raviartthomas12d_HEADERS       raviartthomas12dlocalbasis.hh
                                raviartthomas12dlocalcoefficients.hh
                                raviartthomas12dlocalinterpolation.hh)

# include not needed for CMake
# include $(top_srcdir)/am/global-rules

install(FILES ${raviartthomas12d_HEADERS} DESTINATION ${raviartthomas12ddir})
