// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <iostream>

#include <dune/grid/genericgeometry/referenceelements.hh>

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

  typedef ReferenceElement< Topology, double > RefElement;

  for( unsigned int i = 0; i < RefElement::numCorners; ++i )
  {
    basis.evaluate( p, RefElement::corner( i ), y );
    std::cout << "x = " << RefElement::corner( i ) << ":" << std::endl;
    for( unsigned int i = 0; i < size; ++i )
      std::cout << "    y[ " << i << " ] = " << y[ i ] << std::endl;
  }
}
