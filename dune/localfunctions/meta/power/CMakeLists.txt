set(mydir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/meta/power)

set(my_HEADERS
        basis.hh
        coefficients.hh
        interpolation.hh)

# include not needed for CMake
# include $(top_srcdir)/am/global-rules

install(FILES ${my_HEADERS} DESTINATION ${mydir})
