# which tests to run
set(TESTS
        lagrangeshapefunctiontest
        monomialshapefunctiontest
        virtualshapefunctiontest
        testgenericfem
        test-lagrange1
        test-lagrange2
        test-lagrange3
        test-orthonormal1
        test-orthonormal2
        test-orthonormal3)

set(noinst_HEADERS
        geometries.hh
        test-fe.hh
        test-localfe.hh)

# define the programs
add_executable(testgenericfem testgenericfem.cc)
add_dune_gmp_flags(testgenericfem)
target_link_libraries(testgenericfem ${DUNE_LIBS})

add_executable(lagrangeshapefunctiontest lagrangeshapefunctiontest.cc)
target_link_libraries(lagrangeshapefunctiontest ${DUNE_LIBS})

add_executable(monomialshapefunctiontest monomialshapefunctiontest.cc)
target_link_libraries(monomialshapefunctiontest ${DUNE_LIBS})

add_executable(virtualshapefunctiontest virtualshapefunctiontest.cc)
target_link_libraries(virtualshapefunctiontest ${DUNE_LIBS})

list(APPEND TESTS test-edges0.5)
add_executable(test-edges0.5 test-edges0.5.cc)
target_link_libraries(test-edges0.5 ${DUNE_LIBS})

list(APPEND TESTS test-localfe)
add_executable(test-localfe test-localfe.cc)
target_link_libraries(test-localfe ${DUNE_LIBS})

list(APPEND TESTS test-monomial)
add_executable(test-monomial test-monomial.cc)
target_link_libraries(test-monomial ${DUNE_LIBS})

list(APPEND TESTS test-pk2d)
add_executable(test-pk2d test-pk2d.cc)
target_link_libraries(test-pk2d ${DUNE_LIBS})

list(APPEND TESTS test-power-monomial)
add_executable(test-power-monomial test-power-monomial.cc)
target_link_libraries(test-power-monomial ${DUNE_LIBS})

list(APPEND TESTS test-q1)
add_executable(test-q1 test-q1.cc)
target_link_libraries(test-q1 ${DUNE_LIBS})

list(APPEND TESTS test-q2)
add_executable(test-q2 test-q2.cc)
target_link_libraries(test-q2 ${DUNE_LIBS})

set(LAGRANGE1_TOPOLOGY  "Pyramid< Point >")
set(LAGRANGE1_SETTINGS  TOPOLOGY=${LAGRANGE1_TOPOLOGY})

add_executable(test-lagrange1 test-lagrange.cc)
set_property(TARGET test-lagrange1 APPEND PROPERTY COMPILE_DEFINITIONS ${LAGRANGE1_SETTINGS})
target_link_libraries(test-lagrange1 dunecommon)
add_dune_gmp_flags(test-lagrange1)

set(LAGRANGE2_TOPOLOGY  "Pyramid< ${LAGRANGE1_TOPOLOGY} >")
set(LAGRANGE2_SETTINGS  TOPOLOGY=${LAGRANGE2_TOPOLOGY})

add_executable(test-lagrange2 test-lagrange.cc)
set_property(TARGET test-lagrange2 APPEND PROPERTY COMPILE_DEFINITIONS ${LAGRANGE2_SETTINGS})
target_link_libraries(test-lagrange2 dunecommon)
add_dune_gmp_flags(test-lagrange2)

set(LAGRANGE3_TOPOLOGY  "Pyramid< ${LAGRANGE2_TOPOLOGY} >")
set(LAGRANGE3_SETTINGS  TOPOLOGY=${LAGRANGE3_TOPOLOGY})

add_executable(test-lagrange3 test-lagrange.cc)
set_property(TARGET test-lagrange3 APPEND PROPERTY COMPILE_DEFINITIONS ${LAGRANGE3_SETTINGS})
target_link_libraries(test-lagrange3 dunecommon)
add_dune_gmp_flags(test-lagrange3)

set(ORTHONORMAL1_TOPOLOGY  "Pyramid< Point >")
set(ORTHONORMAL1_SETTINGS  TOPOLOGY=${ORTHONORMAL1_TOPOLOGY})

add_executable(test-orthonormal1 test-orthonormal.cc)
set_property(TARGET test-orthonormal1 APPEND PROPERTY COMPILE_DEFINITIONS ${ORTHONORMAL1_SETTINGS})
target_link_libraries(test-orthonormal1 dunegeometry dunecommon)
add_dune_gmp_flags(test-orthonormal1)

set(ORTHONORMAL2_TOPOLOGY  "Pyramid< ${ORTHONORMAL1_TOPOLOGY} >")
set(ORTHONORMAL2_SETTINGS  TOPOLOGY=${ORTHONORMAL2_TOPOLOGY})

add_executable(test-orthonormal2 test-orthonormal.cc)
set_property(TARGET test-orthonormal2 APPEND PROPERTY COMPILE_DEFINITIONS ${ORTHONORMAL2_SETTINGS})
target_link_libraries(test-orthonormal2 dunegeometry dunecommon)
add_dune_gmp_flags(test-orthonormal2)

set(ORTHONORMAL3_TOPOLOGY  "Prism< ${ORTHONORMAL2_TOPOLOGY} >")
set(ORTHONORMAL3_SETTINGS  TOPOLOGY=${ORTHONORMAL3_TOPOLOGY})

add_executable(test-orthonormal3 test-orthonormal.cc)
set_property(TARGET test-orthonormal3 APPEND PROPERTY COMPILE_DEFINITIONS ${ORTHONORMAL3_SETTINGS})
target_link_libraries(test-orthonormal3 dunegeometry dunecommon)
add_dune_gmp_flags(test-orthonormal3)

list(APPEND TEST test-raviartthomas1simplex)

set(RT1_TOPOLOGY "Pyramid< Point >")
set(RT1_SETTINGS TOPOLOGY=${RT2_TOPOLOGY})
add_executable(test-raviartthomas1simplex test-raviartthomassimplex.cc)
target_link_libraries(test-raviartthomas1simplex ${DUNE_LIBS})
set_property(TARGET test-raviartthomas1simplex APPEND PROPERTY COMPILE_DEFINITIONS
  TOPOLOGY=Pyramid<Point>)
add_dune_gmp_flags(test-raviartthomas1simplex)

list(APPEND TEST test-raviartthomas2simplex)
set(RT2_TOPOLOGY "Pyramid< ${RT1_TOPOLOGY} >")
set(RT2_SETTINGS TOPOLOGY=${RT2_TOPOLOGY})
add_executable(test-raviartthomas2simplex test-raviartthomassimplex.cc)
target_link_libraries(test-raviartthomas2simplex ${DUNE_LIBS})
set_property(TARGET test-raviartthomas2simplex APPEND PROPERTY COMPILE_DEFINITIONS
  ${RT2_SETTINGS})
add_dune_gmp_flags(test-raviartthomas2simplex)

list(APPEND TEST test-raviartthomas3simplex)
set(RT3_TOPOLOGY "Pyramid< ${RT2_TOPOLOGY} >")
set(RT3_SETTINGS TOPOLOGY=${RT3_TOPOLOGY})
add_executable(test-raviartthomas3simplex test-raviartthomassimplex.cc)
target_link_libraries(test-raviartthomas3simplex ${DUNE_LIBS})
set_property(TARGET test-raviartthomas3simplex APPEND PROPERTY COMPILE_DEFINITIONS
  ${RT3_SETTINGS})
add_dune_gmp_flags(test-raviartthomas3simplex)

add_executable(test-raviartthomassimplex-all test-raviartthomassimplex.cc)
target_link_libraries(test-raviartthomassimplex-all ${DUNE_LIBS})
add_dune_gmp_flags(test-raviartthomassimplex-all)

set(TOPOLOGY  "Pyramid< Pyramid< Pyramid< Point > > >")
set(SETTINGS  TOPOLOGY=${TOPOLOGY})

add_executable(test-lagrange-all test-lagrange.cc)
add_dune_gmp_flags(test-lagrange-all)
target_link_libraries(test-lagrange-all dunecommon)

add_executable(test-orthonormal-all test-orthonormal.cc)
add_dune_gmp_flags(test-orthonormal-all)
target_link_libraries(test-orthonormal-all dunecommon)

foreach(i ${TESTS})
  add_test(${i} ${i})
endforeach(i ${TESTS})

# We do not want want to build the tests during make all
# but just build them on demand
add_directory_test_target(_test_target)
add_dependencies(${_test_target} ${TESTS} ${COMPILE_XFAIL_TESTS})
