set(raviartthomas0q2ddir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/raviartthomas/raviartthomas0q2d)

set(raviartthomas0q2d_HEADERS  raviartthomas0q2dall.hh)

# include not needed for CMake
# include $(top_srcdir)/am/global-rules

install(FILES ${raviartthomas0q2d_HEADERS} DESTINATION ${raviartthomas0q2ddir})
