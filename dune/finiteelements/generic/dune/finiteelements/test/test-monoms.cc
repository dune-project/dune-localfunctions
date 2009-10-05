// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <iostream>

#include <dune/finiteelements/monomialbasis.hh>
#include <dune/finiteelements/multiindex.hh>
#include <dune/finiteelements/lagrangepoints.hh>
#include <dune/alglib/multiprecision.hh>

#include <dune/finiteelements/quadrature/genericquadrature.hh>

#ifndef TOPOLOGY
#error "TOPOLOGY not defined."
#endif

using namespace Dune::GenericGeometry;
using namespace Dune;

int main ( int argc, char **argv )
{
  if( argc < 3 )
  {
    std::cerr << "Usage: " << argv[ 0 ] << " <p> <deriv>" << std::endl;
    return 1;
  }

  int p = atoi( argv[ 1 ] );
  int deriv = atoi( argv[ 2 ] );

  typedef TOPOLOGY Topology;
  typedef double Field;

  MonomialBasis< Topology, Field > basis( p );
  const unsigned int derivSize = basis.derivSize( deriv );
  std::vector< Field > y( derivSize * basis.size() );

  typedef Dune::LagrangePoints< Topology, Field > LagrangePoints;
  LagrangePoints points( p );

  std::cout << "Number of base functions:  " << basis.size() << std::endl;
  std::cout << "Number of Lagrange points: " << points.size() << std::endl;

  const LagrangePoints::iterator end = points.end();
  for( LagrangePoints::iterator it = points.begin(); it != end; ++it )
  {
    for( unsigned int i = 0; i < y.size(); ++i )
      y[ i ] = -42.3456789;
    basis.evaluate( deriv, it->point(), &(y[ 0 ]) );
    std::cout << "x = " << field_cast< double >( it->point() )
              << " (codim = " << it->localKey().codim() << ", "
              << "subentity = " << it->localKey().subEntity() << ", "
              << "index = " << it->localKey().index() << "):" << std::endl;
    for( unsigned int i = 0; i < basis.size(); ++i )
    {
      std::cout << "    y[ " << i << " ] = " << y[ i*derivSize ];
      for( unsigned int d = 1; d < derivSize; ++d )
        std::cout << ", " << y[ i*derivSize + d ];
      std::cout << std::endl;
    }
  }
}
