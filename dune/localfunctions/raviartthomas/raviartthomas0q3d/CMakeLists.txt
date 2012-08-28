set(raviartthomas0q3ddir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/raviartthomas/raviartthomas0q3d)

set(raviartthomas0q3d_HEADERS  raviartthomas0q3dall.hh)

# include not needed for CMake
# include $(top_srcdir)/am/global-rules

install(FILES ${raviartthomas0q3d_HEADERS} DESTINATION ${raviartthomas0q3ddir})
