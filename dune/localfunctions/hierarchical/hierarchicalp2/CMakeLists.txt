set(hierarchicalp2dir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/hierarchical/hierarchicalp2)

set(hierarchicalp2_HEADERS  hierarchicalsimplexp2localbasis.hh hierarchicalsimplexp2localinterpolation.hh)

# include not needed for CMake
# include $(top_srcdir)/am/global-rules

install(FILES ${hierarchicalp2_HEADERS} DESTINATION ${hierarchicalp2dir})
