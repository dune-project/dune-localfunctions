set(raviartthomasdir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/raviartthomas)

set(raviartthomas_HEADERS  raviartthomas02d.hh raviartthomas0q2d.hh raviartthomas0q3d.hh raviartthomas0q.hh
                        raviartthomas12d.hh raviartthomas1q2d.hh raviartthomas1q3d.hh raviartthomas2q2d.hh
			raviartthomas0cube2d.hh raviartthomas0cube3d.hh
			raviartthomas1cube2d.hh raviartthomas1cube3d.hh
			raviartthomas2cube2d.hh
			raviartthomassimplex.hh)

set(SUBDIRS
	raviartthomas02d
	raviartthomas0cube2d
	raviartthomas0cube3d
	raviartthomas12d
	raviartthomas1cube2d
	raviartthomas1cube3d
	raviartthomas2cube2d
	raviartthomassimplex)

install(FILES ${raviartthomas_HEADERS} DESTINATION ${raviartthomasdir})
foreach(i ${SUBDIRS})
  if(${i} STREQUAL "test")
    set(opt "EXCLUDE_FROM_ALL")
  endif(${i} STREQUAL "test")
  add_subdirectory(${i} ${opt})
  unset(opt)
endforeach(i ${SUBDIRS})
