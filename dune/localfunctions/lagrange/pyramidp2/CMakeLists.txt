set(pyramidp2dir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/lagrange/pyramidp2)

set(pyramidp2_HEADERS  pyramidp2localbasis.hh pyramidp2localcoefficients.hh pyramidp2localinterpolation.hh)

# include not needed for CMake
# include $(top_srcdir)/am/global-rules

install(FILES ${pyramidp2_HEADERS} DESTINATION ${pyramidp2dir})
