// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <iostream>

#include <dune/finiteelements/lagrangepoints.hh>
#include <dune/finiteelements/monomialbasis.hh>

#ifndef TOPOLOGY
#error "TOPOLOGY not defined."
#endif

using namespace Dune::GenericGeometry;

int main ( int argc, char **argv )
{
  if( argc < 2 )
  {
    std::cerr << "Usage: " << argv[ 0 ] << " <p>" << std::endl;
    return 1;
  }

  int p = atoi( argv[ 1 ] );

  typedef TOPOLOGY Topology;

  Dune::MonomialBasis< Topology, double > basis;
  const unsigned int size = basis.sizes( p )[ p ];
  std::vector< Dune::FieldVector< double, 1 > > y( size );

  typedef Dune::LagrangePoints< Topology, double > LagrangePoints;
  LagrangePoints points( p );

  const LagrangePoints::iterator end = points.end();
  for( LagrangePoints::iterator it = points.begin(); it != end; ++it )
  {
    basis.evaluate( p, *it, y );
    std::cout << "x = " << *it << ":" << std::endl;
    for( unsigned int i = 0; i < size; ++i )
      std::cout << "    y[ " << i << " ] = " << y[ i ] << std::endl;
  }
}
