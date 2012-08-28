set(brezzidouglasmarinidir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/brezzidouglasmarini)

set(brezzidouglasmarini_HEADERS  brezzidouglasmarini12d.hh brezzidouglasmarini1q2d.hh)

set(SUBDIRS  brezzidouglasmarini12d brezzidouglasmarini1q2d)

# include not needed for CMake
# include $(top_srcdir)/am/global-rules

install(FILES ${brezzidouglasmarini_HEADERS} DESTINATION ${brezzidouglasmarinidir})
foreach(i ${SUBDIRS})
  if(${i} STREQUAL "test")
    set(opt "EXCLUDE_FROM_ALL")
  endif(${i} STREQUAL "test")
  add_subdirectory(${i} ${opt})
endforeach(i ${SUBDIRS})
