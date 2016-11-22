dune_add_test(SOURCES testgenericfem.cc)

dune_add_test(SOURCES lagrangeshapefunctiontest.cc)

dune_add_test(SOURCES monomialshapefunctiontest.cc)

dune_add_test(SOURCES virtualshapefunctiontest.cc)

dune_add_test(SOURCES test-edges0.5.cc)

dune_add_test(SOURCES test-localfe.cc)

dune_add_test(SOURCES test-monomial)

dune_add_test(SOURCES test-pk2d.cc)

dune_add_test(SOURCES test-power-monomial.cc)

dune_add_test(SOURCES test-q1.cc)

dune_add_test(SOURCES test-q2.cc)

dune_add_test(NAME test-lagrange1
              SOURCES test-lagrange.cc
              COMPILE_DEFINITIONS TOPOLOGY=Pyramid<Point>)

dune_add_test(NAME test-lagrange2
              SOURCES test-lagrange.cc
              COMPILE_DEFINITIONS "TOPOLOGY=Pyramid<Pyramid<Point> >")

dune_add_test(NAME test-lagrange3
              SOURCES test-lagrange.cc
              COMPILE_DEFINITIONS "TOPOLOGY=Pyramid<Pyramid<Pyramid<Point> > >")

dune_add_test(NAME test-orthonormal1
              SOURCES test-orthonormal.cc
              COMPILE_DEFINITIONS TOPOLOGY=Pyramid<Point>)

dune_add_test(NAME test-orthonormal2
              SOURCES test-orthonormal.cc
              COMPILE_DEFINITIONS "TOPOLOGY=Pyramid<Pyramid<Point> >")

dune_add_test(NAME test-orthonormal3
              SOURCES test-orthonormal.cc
              COMPILE_DEFINITIONS "TOPOLOGY=Pyramid<Pyramid<Pyramid<Point> > >")

dune_add_test(NAME test-raviartthomassimplex1
              SOURCES test-raviartthomassimplex.cc
              COMPILE_DEFINITIONS TOPOLOGY=Pyramid<Point>)

dune_add_test(NAME test-raviartthomassimplex2
              SOURCES test-raviartthomassimplex.cc
              COMPILE_DEFINITIONS "TOPOLOGY=Pyramid<Pyramid<Point> >")

dune_add_test(NAME test-raviartthomassimplex3
              SOURCES test-raviartthomassimplex.cc
              COMPILE_DEFINITIONS "TOPOLOGY=Pyramid<Pyramid<Pyramid<Point> > >")

dune_add_test(NAME test-raviartthomassimplex-all
              SOURCES test-raviartthomassimplex.cc)

dune_add_test(NAME test-lagrange-all
              SOURCES test-lagrange.cc
              TIMEOUT 600)

dune_add_test(NAME test-orthonormal-all
              SOURCES test-orthonormal.cc)
