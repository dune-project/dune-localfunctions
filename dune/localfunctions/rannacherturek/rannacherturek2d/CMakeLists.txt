set(rannacherturek2ddir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/rannacherturek/rannacherturek2d)

set(rannacherturek2d_HEADERS
        rannacherturek2dlocalbasis.hh)

install(FILES ${rannacherturek2d_HEADERS} DESTINATION ${rannacherturek2ddir})
