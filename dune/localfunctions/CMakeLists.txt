set(localfunctionsdir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions)

set(localfunctions_HEADERS
        dualmortarbasis.hh
        lagrange.hh
        mimetic.hh
        monom.hh
        monomial.hh
        orthonormal.hh
        rannacherturek.hh
        raviartthomas.hh
        refined.hh)

set(SUBDIRS
        brezzidouglasmarini
        common
        dualmortarbasis
        hierarchical
        lagrange
        meta
        mimetic
        monomial
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
  unset(opt)
endforeach(i ${SUBDIRS})
