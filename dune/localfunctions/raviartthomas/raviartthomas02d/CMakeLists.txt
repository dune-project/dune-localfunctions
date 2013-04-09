set(raviartthomas02ddir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/raviartthomas/raviartthomas02d)

set(raviartthomas02d_HEADERS  raviartthomas02dlocalbasis.hh raviartthomas02dlocalcoefficients.hh raviartthomas02dlocalinterpolation.hh)

# include not needed for CMake
# include $(top_srcdir)/am/global-rules

install(FILES ${raviartthomas02d_HEADERS} DESTINATION ${raviartthomas02ddir})
