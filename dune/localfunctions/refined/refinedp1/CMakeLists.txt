set(refinedp1dir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/refined/refinedp1)

set(refinedp1_HEADERS  refinedp1localbasis.hh)

# include not needed for CMake
# include $(top_srcdir)/am/global-rules

install(FILES ${refinedp1_HEADERS} DESTINATION ${refinedp1dir})
