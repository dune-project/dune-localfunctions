set(raviartthomas_HEADERS
  raviartthomas02d.hh
  raviartthomas0cube2d.hh
  raviartthomas0cube3d.hh
  raviartthomas12d.hh
  raviartthomas1cube2d.hh
  raviartthomas1cube3d.hh
  raviartthomas2cube2d.hh
  raviartthomas2cube3d.hh
  raviartthomas3cube2d.hh
  raviartthomas4cube2d.hh
  raviartthomascube.hh
  raviartthomassimplex.hh)

set(SUBDIRS
  raviartthomas02d
  raviartthomas0cube2d
  raviartthomas0cube3d
  raviartthomas12d
  raviartthomas1cube2d
  raviartthomas1cube3d
  raviartthomas2cube2d
  raviartthomas2cube3d
  raviartthomas3cube2d
  raviartthomas4cube2d
  raviartthomassimplex)

set(raviartthomasdir ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/raviartthomas)
install(FILES ${raviartthomas_HEADERS} DESTINATION ${raviartthomasdir})

foreach(i ${SUBDIRS})
  if(${i} STREQUAL "test")
    set(opt "EXCLUDE_FROM_ALL")
  endif(${i} STREQUAL "test")
  add_subdirectory(${i} ${opt})
  unset(opt)
endforeach(i ${SUBDIRS})
