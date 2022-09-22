// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#include <config.h>
#include <dune/localfunctions/raviartthomas/raviartthomassimplex/raviartthomassimplexbasis.hh>
#include <dune/localfunctions/utility/field.hh>
#include <dune/localfunctions/utility/basisprint.hh>

/**
 * \file
 * \brief Performs some tests for the generic Raviart-Thomas
 *        shape functions on simplices.
 *
 * The topology can be chosen at compile time by setting TOPOLOGY
 * to a Dune::GeometryType like
 * \code
 * GeometryTypes::simplex(2)
 * \endcode
 * which generates a 2d simplex. If TOPOLOGY is not set, all
 * topologies up to 4d are tested. Note, this may lead to prolonged
 * compiler runs.
 *
 * For debugging purpuse the functions and the derivatives can be
 * printed. You have to define the macro TEST_OUTPUT_FUNCTIONS to
 * activate this function.
 */

#if HAVE_GMP
typedef Dune::GMPField< 128 > StorageField;
typedef Dune::GMPField< 512 > ComputeField;
#else
typedef double StorageField;
typedef double ComputeField;
#endif

template< Dune::GeometryType::Id geometryId >
bool test(unsigned int order)
{
  bool ret = true;

  constexpr Dune::GeometryType geometry = geometryId;
  for (unsigned int o = 0; o <= order; ++o)
  {
    std::cout << "Testing " << geometry << " with order " << o << std::endl;
    typedef Dune::RaviartThomasBasisFactory<geometry.dim(),StorageField,ComputeField> BasisFactory;
    const typename BasisFactory::Object &basis = *BasisFactory::template create<geometry>(o);

    // define the macro TEST_OUTPUT_FUNCTIONS to output files containing functions and
    // derivatives in a human readabible form (aka LaTeX source)
#ifdef TEST_OUTPUT_FUNCTIONS
    std::stringstream name;
    name << "rt_" << geometry << "_p" << o << ".basis";
    std::ofstream out(name.str().c_str());
    Dune::basisPrint<0,BasisFactory,typename BasisFactory::StorageField,geometry>(out,basis);
    Dune::basisPrint<1,BasisFactory,typename BasisFactory::StorageField,geometry>(out,basis);
#endif // TEST_OUTPUT_FUNCTIONS

    // test interpolation
    typedef Dune::RaviartThomasL2InterpolationFactory<geometry.dim(),StorageField> InterpolationFactory;
    const typename InterpolationFactory::Object &interpol = *InterpolationFactory::template create<geometry>(o);
    Dune::LFEMatrix<StorageField> matrix;
    interpol.interpolate(basis,matrix);
    for (unsigned int i=0; i<matrix.rows(); ++i)
      matrix(i,i)-=1;
    for (unsigned int i=0; i<matrix.rows(); ++i)
      for (unsigned int j=0; j<matrix.cols(); ++j)
        if ( std::abs( matrix(i,j) ) > 1000.*Dune::Zero<double>::epsilon() )
          std::cout << "  non-zero entry in interpolation matrix: "
                    << "(" << i << "," << j << ") = " << Dune::field_cast<double>(matrix(i,j))
                    << std::endl;

    BasisFactory::release(&basis);
  }
  if (!ret) {
    std::cout << "   FAILED !" << std::endl;
  }
  return ret;
}

#ifdef CHECKDIM
  #if CHECKDIM==1
      #define CHECKDIM1
  #elif CHECKDIM==2
      #define CHECKDIM2
  #elif CHECKDIM==3
      #define CHECKDIM3
  #elif CHECKDIM==4
      #define CHECKDIM4
  #endif
#else
  #define CHECKDIM1
  #define CHECKDIM2
  #define CHECKDIM3
  #define CHECKDIM4
#endif



int main ( int argc, char **argv )
{
  using namespace Dune;
  using namespace Impl;

  unsigned int order = (argc < 2) ? 5 : atoi(argv[1]);

  if (argc < 2)
  {
    std::cerr << "Usage: " << argv[ 0 ] << " <p>" << std::endl
              << "Using default order of " << order << std::endl;
  }
#ifdef TOPOLOGY
  return (test<TOPOLOGY>(order) ? 0 : 1 );
#else
  bool tests = true;

#ifdef CHECKDIM1
  tests &= test<GeometryTypes::simplex(1)> (order);
#endif

#ifdef CHECKDIM2
  tests &= test<GeometryTypes::simplex(2)> (order);
#endif

#ifdef CHECKDIM3
  tests &= test<GeometryTypes::simplex(3)> (order);
#endif

  // reduce tested order to 4 in 4d unless explicitly asked for more
  if (argc < 2)
    order = 4;

#ifdef CHECKDIM4
  tests &= test<GeometryTypes::simplex(4)> (order);
#endif

  return (tests ? 0 : 1);
#endif // TOPOLOGY
}
