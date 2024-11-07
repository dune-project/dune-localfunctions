# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

dune_add_test(SOURCES bdfmelementtest.cc LINK_LIBRARIES Dune::LocalFunctions)

dune_add_test(SOURCES brezzidouglasmarinielementtest.cc LINK_LIBRARIES Dune::LocalFunctions)

dune_add_test(SOURCES crouzeixraviartelementtest.cc LINK_LIBRARIES Dune::LocalFunctions)

dune_add_test(SOURCES dualmortarelementtest.cc LINK_LIBRARIES Dune::LocalFunctions)

dune_add_test(SOURCES hierarchicalelementtest.cc LINK_LIBRARIES Dune::LocalFunctions)

dune_add_test(SOURCES testgenericfem.cc LINK_LIBRARIES Dune::LocalFunctions)

dune_add_test(SOURCES lagrangeshapefunctiontest.cc LINK_LIBRARIES Dune::LocalFunctions)

dune_add_test(SOURCES monomialshapefunctiontest.cc LINK_LIBRARIES Dune::LocalFunctions)

dune_add_test(SOURCES nedelec1stkindelementtest.cc LINK_LIBRARIES Dune::LocalFunctions)

dune_add_test(SOURCES rannacherturekelementtest.cc LINK_LIBRARIES Dune::LocalFunctions)

dune_add_test(SOURCES raviartthomaselementtest.cc LINK_LIBRARIES Dune::LocalFunctions)

dune_add_test(SOURCES refinedelementtest.cc LINK_LIBRARIES Dune::LocalFunctions)

dune_add_test(SOURCES virtualshapefunctiontest.cc LINK_LIBRARIES Dune::LocalFunctions)

dune_add_test(SOURCES test-edges0.5.cc LINK_LIBRARIES Dune::LocalFunctions)

dune_add_test(SOURCES test-finiteelementcache.cc LINK_LIBRARIES Dune::LocalFunctions)

dune_add_test(SOURCES globalmonomialfunctionstest.cc LINK_LIBRARIES Dune::LocalFunctions)

dune_add_test(SOURCES test-enriched.cc LINK_LIBRARIES Dune::LocalFunctions)

dune_add_test(SOURCES test-pk2d.cc LINK_LIBRARIES Dune::LocalFunctions)

dune_add_test(SOURCES test-power-monomial.cc LINK_LIBRARIES Dune::LocalFunctions)

dune_add_test(SOURCES test-q1.cc LINK_LIBRARIES Dune::LocalFunctions)

dune_add_test(SOURCES test-q2.cc LINK_LIBRARIES Dune::LocalFunctions)

dune_add_test(NAME test-lagrange1
              SOURCES test-lagrange.cc
              LINK_LIBRARIES Dune::LocalFunctions
              COMPILE_DEFINITIONS "CHECKDIM=1")

dune_add_test(NAME test-lagrange2
              SOURCES test-lagrange.cc
              LINK_LIBRARIES Dune::LocalFunctions
              COMPILE_DEFINITIONS "CHECKDIM=2")

dune_add_test(NAME test-lagrange3
              SOURCES test-lagrange.cc
              LINK_LIBRARIES Dune::LocalFunctions
              COMPILE_DEFINITIONS "CHECKDIM=3")

dune_add_test(NAME test-lagrange4
              SOURCES test-lagrange.cc
              LINK_LIBRARIES Dune::LocalFunctions
              COMPILE_DEFINITIONS "CHECKDIM=4"
              TIMEOUT 600)

dune_add_test(NAME test-orthonormal1
              SOURCES test-orthonormal.cc
              LINK_LIBRARIES Dune::LocalFunctions
              COMPILE_DEFINITIONS "CHECKDIM=1")

dune_add_test(NAME test-orthonormal2
              SOURCES test-orthonormal.cc
              LINK_LIBRARIES Dune::LocalFunctions
              COMPILE_DEFINITIONS "CHECKDIM=2")

dune_add_test(NAME test-orthonormal3
              SOURCES test-orthonormal.cc
              LINK_LIBRARIES Dune::LocalFunctions
              COMPILE_DEFINITIONS "CHECKDIM=3")

dune_add_test(NAME test-orthonormal4
              SOURCES test-orthonormal.cc
              LINK_LIBRARIES Dune::LocalFunctions
              COMPILE_DEFINITIONS "CHECKDIM=4")

dune_add_test(NAME test-raviartthomassimplex1
              SOURCES test-raviartthomassimplex.cc
              LINK_LIBRARIES Dune::LocalFunctions
              COMPILE_DEFINITIONS "CHECKDIM=1")

dune_add_test(NAME test-raviartthomassimplex2
              SOURCES test-raviartthomassimplex.cc
              LINK_LIBRARIES Dune::LocalFunctions
              COMPILE_DEFINITIONS "CHECKDIM=2")

dune_add_test(NAME test-raviartthomassimplex3
              SOURCES test-raviartthomassimplex.cc
              LINK_LIBRARIES Dune::LocalFunctions
              COMPILE_DEFINITIONS "CHECKDIM=3")

dune_add_test(NAME test-raviartthomassimplex4
              SOURCES test-raviartthomassimplex.cc
              LINK_LIBRARIES Dune::LocalFunctions
              COMPILE_DEFINITIONS "CHECKDIM=4")

dune_add_test(NAME test-nedelecsimplex2
              SOURCES test-nedelecsimplex.cc
              LINK_LIBRARIES Dune::LocalFunctions
              COMPILE_DEFINITIONS "CHECKDIM=2")

dune_add_test(NAME test-nedelecsimplex3
              SOURCES test-nedelecsimplex.cc
              LINK_LIBRARIES Dune::LocalFunctions
              COMPILE_DEFINITIONS "CHECKDIM=3")
