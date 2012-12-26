# $Id: Makefile.am 4150 2008-05-15 16:12:46Z christi $

if( dune_grid_FOUND )
if( ALUGrid_FOUND )
  set(ALUPROG  refinedp1test)
endif( ALUGrid_FOUND )
endif( dune_grid_FOUND )

# which tests to run
set(TESTS  lagrangeshapefunctiontest
        monomshapefunctiontest
        virtualshapefunctiontest
        testgenericfem
        ${ALUPROG})


# programs just to build when "make check" is used
set(check_PROGRAMS  ${TESTS})

# output code coverage
#AM_CXXFLAGS = @AM_CPPFLAGS@ -fprofile-arcs -ftest-coverage

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
add_executable(monomshapefunctiontest monomshapefunctiontest.cc)
target_link_libraries(monomshapefunctiontest ${DUNE_LIBS})

add_executable(virtualshapefunctiontest virtualshapefunctiontest.cc)
target_link_libraries(virtualshapefunctiontest ${DUNE_LIBS})

if( dune-grid_FOUND)
if( ALUGrid_FOUND)
add_executable(refinedp1test refinedp1test.cc)
dune_add_alugrid_flags(efinedp1test)
endif( ALUGrid_FOUND )
endif( dune-grid_FOUND )

list(APPEND TESTS test-edges0.5)
add_executable(test-edges0.5 test-edges0.5.cc)
target_link_libraries(test-edges0.5 ${DUNE_LIBS})

list(APPEND TESTS test-localfe)
add_executable(test-localfe test-localfe.cc)
target_link_libraries(test-localfe ${DUNE_LIBS})

list(APPEND TESTS test-monom)
add_executable(test-monom test-monom.cc)
target_link_libraries(test-monom ${DUNE_LIBS})

list(APPEND TESTS test-pk2d)
add_executable(test-pk2d test-pk2d.cc)
target_link_libraries(test-pk2d ${DUNE_LIBS})

list(APPEND TESTS test-power-monom)
add_executable(test-power-monom test-power-monom.cc)
target_link_libraries(test-power-monom ${DUNE_LIBS})

list(APPEND TESTS test-q1)
add_executable(test-q1 test-q1.cc)
target_link_libraries(test-q1 ${DUNE_LIBS})

list(APPEND TESTS test-q2)
add_executable(test-q2 test-q2.cc)
target_link_libraries(test-q2 ${DUNE_LIBS})

list(APPEND TEST test-raviartthomas1simplex)
add_executable(test-raviartthomas1simplex test-raviartthomassimplex.cc)
target_link_libraries(test-raviartthomas1simplex ${DUNE_LIBS})
set_property(TARGET test-raviartthomas1simplex APPEND PROPERTY COMPILE_DEFINITIONS
  TOPOLOGY=Pyramid<Point>)
add_dune_gmp_flags(test-raviartthomas1simplex)

list(APPEND TEST test-raviartthomas2simplex)
add_executable(test-raviartthomas2simplex test-raviartthomassimplex.cc)
target_link_libraries(test-raviartthomas2simplex ${DUNE_LIBS})
set_property(TARGET test-raviartthomas2simplex APPEND PROPERTY COMPILE_DEFINITIONS
  "TOPOLOGY=\"Pyramid<Pyramid<Point> >\"")
add_dune_gmp_flags(test-raviartthomas2simplex)

list(APPEND TEST test-raviartthomas3simplex)
add_executable(test-raviartthomas3simplex test-raviartthomassimplex.cc)
target_link_libraries(test-raviartthomas1simplex ${DUNE_LIBS})
set_property(TARGET test-raviartthomas1simplex APPEND PROPERTY COMPILE_DEFINITIONS
  "TOPOLOGY=\"Pyramid<Pyramid<Pyramid<Point> > >\"")
add_dune_gmp_flags(test-raviartthomas3simplex)

add_executable(test-raviartthomassimplex-all test-raviartthomassimplex.cc)
target_link_libraries(test-raviartthomassimplex-all ${DUNE_LIBS})
add_dune_gmp_flags(test-raviartthomassimplex-all)

foreach(i ${TESTS})
  add_test(${i} ${i})
endforeach(i ${TESTS})

# We do not want want to build the tests during make all
# but just build them on demand
add_directory_test_target(_test_target)
add_dependencies(${_test_target} ${TESTS} ${COMPILE_XFAIL_TESTS})
