set(rannacherturek3ddir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/rannacherturek/rannacherturek3d)

set(rannacherturek3d_HEADERS
        rannacherturek3dlocalbasis.hh)

# include not needed for CMake
# include $(top_srcdir)/am/global-rules

install(FILES ${rannacherturek3d_HEADERS} DESTINATION ${rannacherturek3ddir})
