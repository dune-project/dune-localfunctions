set(SUBDIRS  p0 p1 p23d pk1d pk2d pk3d prismp1 prismp2 pyramidp1 pyramidp2 q1 q2 qk)

set(lagrangedir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/lagrange)
set(lagrange_HEADERS  emptypoints.hh        lagrangebasis.hh
                   equidistantpoints.hh  lagrangecoefficients.hh
                   interpolation.hh
                   p0.hh p1.hh p23d.hh p2.hh pk1d.hh pk2d.hh pk3d.hh pk.hh pq22d.hh
                   prismp1.hh prismp2.hh pyramidp1.hh pyramidp2.hh q1.hh q2.hh
                   qk.hh pqkfactory.hh)

set(check_PROGRAMS  test-lagrange-all)
set(EXTRA_PROGRAMS  test-lagrange)


set(TOPOLOGY  "Pyramid< Pyramid< Pyramid< Point > > >")
set(SETTINGS  TOPOLOGY="${TOPOLOGY}")

add_executable(test-lagrange EXCLUDE_FROM_ALL test.cc)
set_property(TARGET test-lagrange APPEND PROPERTY COMPILE_DEFINTIONS ${SETTING})
target_link_libraries(test-lagrange dunecommon)
dune_add_gmp_flags(test-lagrange)
add_executable(test-lagrange-all EXCLUDE_FROM_ALL test.cc)
set_property(TARGET test-lagrange-all APPEND PROPERTY COMPILE_DEFINTIONS ${SETTING})
dune_add_gmp_flags(test-lagrange-all)
target_link_libraries(test-lagrange-all dunecommon)

install(FILES ${lagrange_HEADERS} DESTINATION ${lagrangedir})
foreach(i ${SUBDIRS})
  if(${i} STREQUAL "test")
    set(opt "EXCLUDE_FROM_ALL")
  endif(${i} STREQUAL "test")
  add_subdirectory(${i} ${opt})
endforeach(i ${SUBDIRS})
