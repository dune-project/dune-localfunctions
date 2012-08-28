set(rannacherturekdir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/rannacherturek)

set(rannacherturek_HEADERS  rannacherturek2d.hh)

set(SUBDIRS  rannacherturek2d)

# include not needed for CMake
# include $(top_srcdir)/am/global-rules

install(FILES ${rannacherturek_HEADERS} DESTINATION ${rannacherturekdir})
foreach(i ${SUBDIRS})
  if(${i} STREQUAL "test")
    set(opt "EXCLUDE_FROM_ALL")
  endif(${i} STREQUAL "test")
  add_subdirectory(${i} ${opt})
endforeach(i ${SUBDIRS})
