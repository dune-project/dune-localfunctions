set(pyramidp1dir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/lagrange/pyramidp1)

set(pyramidp1_HEADERS  pyramidp1localbasis.hh pyramidp1localcoefficients.hh pyramidp1localinterpolation.hh)

# include not needed for CMake
# include $(top_srcdir)/am/global-rules

install(FILES ${pyramidp1_HEADERS} DESTINATION ${pyramidp1dir})
