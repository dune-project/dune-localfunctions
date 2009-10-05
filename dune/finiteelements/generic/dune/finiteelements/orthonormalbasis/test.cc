// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <dune/finiteelements/orthonormalbasis.hh>

#include <dune/grid/quadrature/genericquadrature.hh>

using namespace Dune;
using namespace GenericGeometry;

int main ( int argc, char **argv )
{
  if( argc < 2 )
  {
    std::cerr << "Usage: " << argv[ 0 ] << " <p>" << std::endl;
    return 1;
  }

  typedef TOPOLOGY Topology;

  const unsigned int order = atoi( argv[ 1 ] );
  OrthonormalBasis< Topology, double > basis( order );

  const unsigned int size = basis.size( order );
  std::vector< FieldVector< double, 1 > > y( size );

  std::vector< FieldVector< double, 1 > > m( size * size );
  for( unsigned int i = 0; i < size * size; ++i )
    m[ i ] = 0;

  typedef QuadratureRuleImpl< Topology > Quadrature;
  Quadrature quadrature( 2*order+1 );
  const Quadrature::iterator end = quadrature.end();
  for( Quadrature::iterator it = quadrature.begin(); it != end; ++it )
  {
    basis.evaluate( order, it->position(), y );
    std::cout << "Point: " << it->position() << ", " << std::endl;
    std::cout << "weight: " << it->weight() << ":" << std::endl;
    for( unsigned int i = 0; i < size; ++i )
    {
      std::cout << "-> " << y[ i ] << std::endl;
      for( unsigned int j = 0; j < size; ++j )
        m[ i*size + j ] += it->weight() * y[ i ] * y[ j ];
    }
  }

  for( unsigned int i = 0; i < size; ++i )
  {
    for( unsigned int j = 0; j < size; ++j )
    {
      std::cout << "i = " << i << ", j = " << j << ": " << m[ i*size + j ] << std::endl;
    }
  }
}
