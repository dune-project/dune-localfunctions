// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#include <config.h>

#include <dune/geometry/type.hh>

#include <dune/localfunctions/utility/field.hh>
#include <dune/localfunctions/utility/basisprint.hh>

#include <dune/localfunctions/lagrange/equidistantpoints.hh>
#include <dune/localfunctions/lagrange/lagrangebasis.hh>

/**
 * \file
 * \brief Performs some tests for the generic Lagrange
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

template <class Basis,class Points>
bool test(const Basis &basis, const Points &points, bool verbose)
{
  bool ret = true;
  std::vector< Dune::FieldVector< double, 1 > > y( basis.size() );
  for( unsigned int index = 0; index < points.size(); ++index )
  {
    if (verbose)
      std::cout << index << "   " << points[ index ].point() << " "
                << points[ index ].localKey()
                << std::endl;
    basis.evaluate( points[ index ].point(), y );
    bool first = true;
    for( unsigned int i = 0; i < y.size(); ++i )
    {
      if( fabs( y[ i ] - double( i == index ) ) > 1e-10 )
      {
        if (first) {
          std::cout << "ERROR: "
                    << index << " -> "
                    << "x = " << points[ index ].point()
                    << " (codim = " << points[ index ].localKey().codim() << ", "
                    << "subentity = " << points[ index ].localKey().subEntity() << ", "
                    << "index = " << points[ index ].localKey().index() << "):" << std::endl;
          first = false;
        }
        if (1)
          std::cout << "         y[ " << i << " ] = " << y[ i ] << " "
                    << "         error : " << fabs( y[ i ] - double( i == index ) )
                    << std::endl;
        ret = false;
      }
    }
  }
  return ret;
}

template <Dune::GeometryType::Id geometryId>
bool test(unsigned int order, bool verbose = false)
{
  constexpr Dune::GeometryType geometry = geometryId;
  typedef Dune::LagrangeBasisFactory<Dune::EquidistantPointSet, geometry.dim(), StorageField, ComputeField> BasisFactory;
  typedef Dune::LagrangeCoefficientsFactory< Dune::EquidistantPointSet, geometry.dim(), double > LagrangeCoefficientsFactory;

  bool ret = true;

  for (unsigned int o = 0; o <= order; ++o)
  {
    const typename LagrangeCoefficientsFactory::Object *pointsPtr = LagrangeCoefficientsFactory::template create< geometry >( o );

    if ( pointsPtr == 0)
      continue;

    std::cout << "Testing " << geometry << " with order " << o << std::endl;

    typename BasisFactory::Object &basis = *BasisFactory::template create<geometry>(o);

    ret |= test(basis,*pointsPtr,verbose);

    // define the macro TEST_OUTPUT_FUNCTIONS to output files containing functions and
    // derivatives in a human readabible form (aka LaTeX source)
#ifdef TEST_OUTPUT_FUNCTIONS
    std::stringstream name;
    name << "lagrange_" << geometry << "_p" << o << ".basis";
    std::ofstream out(name.str().c_str());
    Dune::basisPrint<0,BasisFactory,typename BasisFactory::StorageField,geometry>(out,basis);
    Dune::basisPrint<1,BasisFactory,typename BasisFactory::StorageField,geometry>(out,basis);
#endif // TEST_OUTPUT_FUNCTIONS

    LagrangeCoefficientsFactory::release( pointsPtr );
    BasisFactory::release( &basis );
  }

  if (verbose)
    std::cout << std::endl << std::endl << std::endl;
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
  tests &= test<GeometryTypes::cube(1)> (order);
  tests &= test<GeometryTypes::simplex(1)> (order);
#endif

#ifdef CHECKDIM2
  tests &= test<GeometryTypes::cube(2)> (order);
  tests &= test<GeometryTypes::simplex(2)> (order);
#endif

#ifdef CHECKDIM3
  tests &= test<GeometryTypes::cube(3)> (order);
  tests &= test<GeometryTypes::prism> (order);
  tests &= test<GeometryTypes::pyramid> (order);
  tests &= test<GeometryTypes::simplex(3)> (order);
#endif

  // reduce tested order to 4 in 4d unless explicitly asked for more
  if (argc < 2)
    order = 4;

#ifdef CHECKDIM4
  tests &= test<GeometryTypes::cube(4)> (order);
  tests &= test<GeometryTypes::simplex(4)> (order);
#endif

  return (tests ? 0 : 1);
#endif // TOPOLOGY
}
