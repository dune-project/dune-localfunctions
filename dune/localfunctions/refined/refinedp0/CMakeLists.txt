set(refinedp0dir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/refined/refinedp0)

set(refinedp0_HEADERS  refinedp0localbasis.hh refinedp0localcoefficients.hh refinedp0localinterpolation.hh)

# include not needed for CMake
# include $(top_srcdir)/am/global-rules

install(FILES ${refinedp0_HEADERS} DESTINATION ${refinedp0dir})
