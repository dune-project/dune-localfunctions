set(orthonormaldir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/orthonormal)
set(orthonormal_HEADERS  orthonormalbasis.hh orthonormalcompute.hh)
set(check_PROGRAMS  test-orthonormal)
set(EXTRA_PROGRAMS  test-orthonormal-all)

set(TOPOLOGY  "Prism< Pyramid< Pyramid< Point > > >")
set(SETTINGS  TOPOLOGY=${TOPOLOGY})

add_executable(test-orthonormal EXCLUDE_FROM_ALL test.cc)
set_property(TARGET test-orthonormal APPEND PROPERTY COMPILE_DEFINITIONS ${SETTINGS})
add_dune_gmp_flags(test-orthonormal)
target_link_libraries(test-orthonormal dunegeometry  dunecommon)
add_executable(test-orthonormal-all EXCLUDE_FROM_ALL test.cc)
set_property(TARGET test-orthonormal-all APPEND PROPERTY COMPILE_DEFINITIONS ${SETTINGS})
add_dune_gmp_flags(test-orthonormal-all)
target_link_libraries(test-orthonormal-all dunegeometry dunecommon)

install(FILES ${orthonormal_HEADERS} DESTINATION ${orthonormaldir})
