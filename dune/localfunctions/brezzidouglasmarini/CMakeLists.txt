set(brezzidouglasmarinidir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/brezzidouglasmarini)

set(brezzidouglasmarini_HEADERS
  brezzidouglasmarini1cube2d.hh
  brezzidouglasmarini2cube2d.hh
  brezzidouglasmarini1q2d.hh
  brezzidouglasmarini1simplex2d.hh
  brezzidouglasmarini12d.hh)

set(SUBDIRS
  brezzidouglasmarini1cube2d
  brezzidouglasmarini2cube2d
  brezzidouglasmarini1simplex2d)

install(FILES ${brezzidouglasmarini_HEADERS} DESTINATION ${brezzidouglasmarinidir})
foreach(i ${SUBDIRS})
  if(${i} STREQUAL "test")
    set(opt "EXCLUDE_FROM_ALL")
  endif(${i} STREQUAL "test")
  add_subdirectory(${i} ${opt})
endforeach(i ${SUBDIRS})
