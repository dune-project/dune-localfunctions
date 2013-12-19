set(SUBDIRS  p0 p1 p23d pk1d pk2d pk3d prismp1 prismp2 pyramidp1 pyramidp2 q1 q2 qk)

set(lagrangedir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/lagrange)
set(lagrange_HEADERS  emptypoints.hh        lagrangebasis.hh
                   equidistantpoints.hh  lagrangecoefficients.hh
                   interpolation.hh
                   p0.hh p1.hh p23d.hh p2.hh pk1d.hh pk2d.hh pk3d.hh pk.hh pq22d.hh
                   prismp1.hh prismp2.hh pyramidp1.hh pyramidp2.hh q1.hh q2.hh
                   qk.hh pqkfactory.hh)

install(FILES ${lagrange_HEADERS} DESTINATION ${lagrangedir})
foreach(i ${SUBDIRS})
  if(${i} STREQUAL "test")
    set(opt "EXCLUDE_FROM_ALL")
  endif(${i} STREQUAL "test")
  add_subdirectory(${i} ${opt})
  unset(opt)
endforeach(i ${SUBDIRS})
