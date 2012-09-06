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
dune_add_gmp_flags(testgenericfem)
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

foreach(i ${TESTS})
  add_test(${i} ${i})
endforeach(i ${TESTS})

# We do not want want to build the tests during make all
# but just build them on demand
add_directory_test_target(_test_target)
add_dependencies(${_test_target} ${TESTS} ${COMPILE_XFAIL_TESTS})
