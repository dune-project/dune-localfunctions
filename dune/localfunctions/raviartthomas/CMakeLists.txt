set(raviartthomasdir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/raviartthomas)

set(raviartthomas_HEADERS  raviartthomas02d.hh raviartthomas0q2d.hh raviartthomas0q3d.hh raviartthomas0q.hh
                        raviartthomas12d.hh raviartthomas1q2d.hh raviartthomas1q3d.hh raviartthomas2q2d.hh
                        interpolation.hh  raviartthomasbasis.hh  raviartthomasprebasis.hh)

set(SUBDIRS  raviartthomas02d raviartthomas0q2d raviartthomas0q3d
          raviartthomas12d raviartthomas1q2d raviartthomas1q3d
          raviartthomas2q2d)

set(check_PROGRAMS  test-raviartthomas)
set(EXTRA_PROGRAMS  test-raviartthomas-all)

set(TOPOLOGY  "Pyramid< Pyramid< Point > >")
set(SETTINGS  TOPOLOGY=${TOPOLOGY})

add_executable(test-raviartthomas test.cc)
set_property(TARGET test-raviartthomas APPEND PROPERTY COMPILE_DEFINITIONS  ${SETTINGS})
dune_add_gmp_flags(test-raviartthomas)
target_link_libraries(test-raviartthomas  dunegeometry dunecommon)
add_executable(test-raviartthomas-all test.cc)
dune_add_gmp_flags(test-raviartthomas-all)
target_link_libraries(test-raviartthomas-all  dunegeometry dunecommon)

# include not needed for CMake
# include $(top_srcdir)/am/global-rules

message(AUTHOR_WARNING "TODO: tests whether check_PROGRAM is treated correctly")
#install(TARGETS ${check_PROGRAMS} DESTINATION ${checkdir})
install(FILES ${raviartthomas_HEADERS} DESTINATION ${raviartthomasdir})
foreach(i ${SUBDIRS})
  if(${i} STREQUAL "test")
    set(opt "EXCLUDE_FROM_ALL")
  endif(${i} STREQUAL "test")
  add_subdirectory(${i} ${opt})
endforeach(i ${SUBDIRS})
