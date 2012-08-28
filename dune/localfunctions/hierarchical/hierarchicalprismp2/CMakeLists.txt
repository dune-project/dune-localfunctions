set(hierarchicalprismp2dir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/hierarchical/hierarchicalprismp2)

set(hierarchicalprismp2_HEADERS  hierarchicalprismp2localbasis.hh hierarchicalprismp2localinterpolation.hh)

# include not needed for CMake
# include $(top_srcdir)/am/global-rules

install(FILES ${hierarchicalprismp2_HEADERS} DESTINATION ${hierarchicalprismp2dir})
