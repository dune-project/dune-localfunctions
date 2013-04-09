set(hierarchicalp2withelementbubbledir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/hierarchical/hierarchicalp2withelementbubble)

set(hierarchicalp2withelementbubble_HEADERS  hierarchicalsimplexp2withelementbubble.hh)

# include not needed for CMake
# include $(top_srcdir)/am/global-rules

install(FILES ${hierarchicalp2withelementbubble_HEADERS} DESTINATION ${hierarchicalp2withelementbubbledir})
