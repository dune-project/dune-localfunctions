set(localfunctionsdir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions)

set(localfunctions_HEADERS
        lagrange.hh
        mimetic.hh
        monom.hh
        orthonormal.hh
        rannacherturek.hh
        raviartthomas.hh
        refined.hh)

set(SUBDIRS
        brezzidouglasmarini
        common
        hierarchical
        lagrange
        meta
        mimetic
        monom
        orthonormal
        rannacherturek
        raviartthomas
        refined
        test
        utility
        whitney)

# include not needed for CMake
# include $(top_srcdir)/am/global-rules

install(FILES ${localfunctions_HEADERS} DESTINATION ${localfunctionsdir})
foreach(i ${SUBDIRS})
  if(${i} STREQUAL "test")
    set(opt "EXCLUDE_FROM_ALL")
  endif(${i} STREQUAL "test")
  add_subdirectory(${i} ${opt})
endforeach(i ${SUBDIRS})
