set(mimeticdir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/mimetic)

set(mimetic_HEADERS  mimeticall.hh)

# include not needed for CMake
# include $(top_srcdir)/am/global-rules

install(FILES ${mimetic_HEADERS} DESTINATION ${mimeticdir})
