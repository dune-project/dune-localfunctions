set(commondir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/refined/common)

set(common_HEADERS  refinedsimplexlocalbasis.hh)

# include not needed for CMake
# include $(top_srcdir)/am/global-rules

install(FILES ${common_HEADERS} DESTINATION ${commondir})
