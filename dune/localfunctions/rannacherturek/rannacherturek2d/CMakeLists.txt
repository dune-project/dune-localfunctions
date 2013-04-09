set(rannacherturek2ddir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/rannacherturek/rannacherturek2d)

set(rannacherturek2d_HEADERS
        rannacherturek2dlocalbasis.hh
        rannacherturek2dlocalcoefficients.hh
        rannacherturek2dlocalinterpolation.hh)

# include not needed for CMake
# include $(top_srcdir)/am/global-rules

install(FILES ${rannacherturek2d_HEADERS} DESTINATION ${rannacherturek2ddir})
