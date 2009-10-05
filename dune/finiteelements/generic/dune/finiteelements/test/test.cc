// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <iostream>
#include <vector>

#include <dune/alglib/multiprecision.hh>
#include <dune/finiteelements/quadrature/genericquadrature.hh>
#include <dune/finiteelements/lagrangepoints.hh>
#include <dune/finiteelements/monomialbasis.hh>
#include <dune/finiteelements/multiindex.hh>

#ifndef TOPOLOGY
#error "TOPOLOGY not defined."
#endif

using namespace Dune;
using namespace GenericGeometry;

template <class Topology>
void lagrangePointTest(unsigned int p)
{
  typedef double Field;

  typedef MonomialBasis< Topology, Field > Basis;
  Basis basis(p);

  unsigned int size = basis.size();
  std::cout << "Number of base functions:  " << size << std::endl;

  typedef LagrangePoints< Topology, Field > LagrangePoints;
  LagrangePoints points( p );
  std::cout << "Number of Lagrange points: " << points.size() << std::endl;
  std::vector< Field > y( size );

  const typename LagrangePoints::iterator end = points.end();
  for( typename LagrangePoints::iterator it = points.begin(); it != end; ++it )
  {
    basis.evaluate( it->point(), &(y[0]) );
    std::cout << "x = " << field_cast<double>(it->point())
              << " (codim = " << it->localKey().codim() << ", "
              << "subentity = " << it->localKey().subEntity() << ", "
              << "index = " << it->localKey().index() << "):" << std::endl;
    for( unsigned int i = 0; i < size; ++i )
      std::cout << "    y[ " << i << " ] = " << field_cast<double>(y[ i ]) << std::endl;
  }
}
template <class Topology>
void quadratureTest(unsigned int p)
{
  typedef double Field;

  typedef MonomialBasis< Topology, Field > Basis;
  Basis basis(p);

  unsigned int size = basis.size();
  std::cout << "Number of base functions:  " << size << std::endl;
  std::cout << std::endl;
  std::cout << ">>> Testing quadrature of order " << (2*p+1) << "..." << std::endl;

  std::vector< FieldVector< Field, 1 > > yquad( size );
  for( unsigned int i = 0; i < size; ++i )
    yquad[ i ] = 0;

  std::vector< Field > y( size );
  GenericQuadrature< Topology,Field > quadrature( 2*p+1 );
  const unsigned int quadratureSize = quadrature.size();
  for( unsigned int qi = 0; qi < quadratureSize; ++qi )
  {
    basis.evaluate( quadrature.point( qi ), y );
    std::cout << "x = " << quadrature.point( qi ) << ":" << std::endl;
    for( unsigned int i = 0; i < size; ++i )
    {
      yquad[ i ] += quadrature.weight( qi ) * y[ i ];
      std::cout << "    y[ " << i << " ] = " << y[ i ] << std::endl;
    }
  }

  std::vector< double > yint( size );
  basis.integral( yint );
  for( unsigned int i = 0; i < size; ++i )
  {
    if( fabs( yquad[ i ] - yint[ i ] ) < 1e-10 )
      continue;
    std::cerr << "Quadrature and Integral differ for basis function " << i << "." << std::endl;
    std::cout << "    quadrature: " << yquad[ i ] << std::endl;
    std::cout << "    integral:   " << yint[ i ] << std::endl;
  }
}
template <class Topology>
void multiIndexTest(unsigned int p)
{
  const int dimension = Topology::dimension;
  typedef MultiIndex< dimension > Field;

  typedef MonomialBasis< Topology, Field > Basis;
  Basis basis(p);

  unsigned int size = basis.size();
  std::cout << "Number of base functions:  " << size << std::endl;

  std::cout << ">>> Polynomial representation of the basis functions:" << std::endl;
  FieldVector< Field, dimension > x;
  for( int i = 0; i < dimension; ++i )
    x[ i ].set( i, 1 );
  std::vector< Field > val( size );
  basis.evaluate( x, val );
  std::cout << val << std::endl;

  std::cout << ">>> integral of basis functions:" << std::endl;
  std::vector< Field > integral( size );
  basis.integral( integral );
  std::cout << integral << std::endl;
}

int main ( int argc, char **argv )
{
  typedef TOPOLOGY Topology;
  if( argc < 2 )
  {
    std::cerr << "Usage: " << argv[ 0 ] << " <p>" << std::endl;
    return 1;
  }
  int p = atoi( argv[ 1 ] );

  lagrangePointTest<Topology>(p);
  quadratureTest<Topology>(p);
  multiIndexTest<Topology>(p);
}
