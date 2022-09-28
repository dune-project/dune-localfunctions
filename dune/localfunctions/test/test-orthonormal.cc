// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#include <config.h>

#include <dune/geometry/type.hh>
#include <dune/geometry/quadraturerules.hh>

#include <dune/localfunctions/utility/field.hh>
#include <dune/localfunctions/utility/basisprint.hh>
#include <dune/localfunctions/orthonormal/orthonormalbasis.hh>

/**
 * \file
 * \brief Performs some tests for the generic orthonormal
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
    typedef Dune::OrthonormalBasisFactory<geometry.dim(),StorageField,ComputeField> BasisFactory;
    const typename BasisFactory::Object &basis = *BasisFactory::template create<geometry>(o);

    const unsigned int size = basis.size( );

    std::vector< Dune::FieldVector< double, 1 > > y( size );

    std::vector< Dune::FieldVector< double, 1 > > m( size * size );
    for( unsigned int i = 0; i < size * size; ++i )
      m[ i ] = 0;

    const Dune::QuadratureRule<double,geometry.dim()> &quadrature =
      Dune::QuadratureRules<double,geometry.dim()>::rule(geometry,2*order+1);
    const unsigned int quadratureSize = quadrature.size();
    for( unsigned int qi = 0; qi < quadratureSize; ++qi )
    {
      basis.evaluate( quadrature[qi].position(), y );
      for( unsigned int i = 0; i < size; ++i )
      {
        for( unsigned int j = 0; j < size; ++j )
          m[ i*size + j ] += quadrature[qi].weight() * y[ i ] * y[ j ];
      }
    }

    for( unsigned int i = 0; i < size; ++i )
    {
      for( unsigned int j = 0; j < size; ++j )
      {
        const double value = m[ i*size + j ];
        if( std::abs( value - double( i == j ) ) > 1200.*Dune::Zero<double>::epsilon() ) {
          std::cout << "i = " << i << ", j = " << j << ": " << std::abs( value - double( i == j ) ) << std::endl;
          ret = false;
        }
      }
    }

    // define the macro TEST_OUTPUT_FUNCTIONS to output files containing functions and
    // derivatives in a human readabible form (aka LaTeX source)
#ifdef TEST_OUTPUT_FUNCTIONS
    std::stringstream name;
    name << "orthonormal_" << geometry << "_p" << o << ".basis";
    std::ofstream out(name.str().c_str());
    Dune::basisPrint<0,BasisFactory,typename BasisFactory::StorageField,geometry>(out,basis);
    Dune::basisPrint<1,BasisFactory,typename BasisFactory::StorageField,geometry>(out,basis);
#endif // TEST_OUTPUT_FUNCTIONS

    BasisFactory::release(&basis);
  }
  if (!ret) {
    std::cout << "   FAILED !" << std::endl;
  }
  std::cout << std::endl;
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

  const unsigned int order = (argc < 2) ? 5 : atoi(argv[1]);

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

#ifdef CHECKDIM4
  tests &= test<GeometryTypes::cube(4)> (order);
  tests &= test<GeometryTypes::simplex(4)> (order);
#endif

  return (tests ? 0 : 1);
#endif // TOPOLOGY
}
