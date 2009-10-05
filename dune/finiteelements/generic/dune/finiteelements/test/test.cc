// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <iostream>

#include <dune/finiteelements/lagrangepoints.hh>
#include <dune/finiteelements/monomialbasis.hh>
#include <dune/finiteelements/multiindex.hh>

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
  const int dimension = Topology::dimension;

  Dune::MonomialBasis< Topology, double > basis;
  const unsigned int size = basis.size( p );
  std::vector< Dune::FieldVector< double, 1 > > y( size );

  typedef Dune::LagrangePoints< Topology, double > LagrangePoints;
  LagrangePoints points( p );

  std::cout << "Number of base functions:  " << size << std::endl;
  std::cout << "Number of Lagrange points: " << points.size() << std::endl;

  const LagrangePoints::iterator end = points.end();
  for( LagrangePoints::iterator it = points.begin(); it != end; ++it )
  {
    basis.evaluate( p, it->point(), y );
    std::cout << "x = " << it->point()
              << " (codim = " << it->localKey().codim() << ", "
              << "subentity = " << it->localKey().subentity() << ", "
              << "index = " << it->localKey().index() << "):" << std::endl;
    for( unsigned int i = 0; i < size; ++i )
      std::cout << "    y[ " << i << " ] = " << y[ i ] << std::endl;
  }

  if( false )
  {
    typedef Dune::StandardBiMonomialBasis< 3,double > Basis;
    Basis basis;
    const unsigned int size = basis.size( p );
    std::vector< Dune::FieldVector< double, 1 > > y( size );

    typedef Dune::LagrangePoints< Basis::Topology, double > LagrangePoints;
    LagrangePoints points( p );

    const LagrangePoints::iterator end = points.end();
    for( LagrangePoints::iterator it = points.begin(); it != end; ++it )
    {
      basis.evaluate( p, it->point(), y );
      std::cout << "x = " << it->point() << ":" << std::endl;
      for( unsigned int i = 0; i < size; ++i )
        std::cout << "    y[ " << i << " ] = " << y[ i ] << std::endl;
    }
  }

  if( true )
  {
    std::cout << std::endl;
    std::cout << "polynomial representation of the basis functions:" << std::endl;
    typedef Dune::MultiIndex< dimension > MultiIndex;
    Dune::MonomialBasis< Topology, MultiIndex  > basis;
    const unsigned int size = basis.size( p );
    std::vector< Dune::FieldVector< MultiIndex, 1 > > y( size );
    Dune::FieldVector< MultiIndex, dimension > x;
    for( int i = 0; i < dimension; ++i )
      x[ i ].set( i, 1 );
    basis.evaluate( p, x, y );
    for( size_t i = 0; i < y.size(); ++i )
      std::cout << y[ i ] << std::endl;
  }
}
