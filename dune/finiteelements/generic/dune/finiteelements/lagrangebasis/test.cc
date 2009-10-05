// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <dune/finiteelements/lagrangebasis.hh>
#include <dune/finiteelements/quadrature/genericquadrature.hh>

using namespace Dune;
using namespace GenericGeometry;

template <class Topology>
bool test(unsigned int order) {
  typedef AlgLib::MultiPrecision<128> StorageField;
  std::cout << "Testing " << Topology::name() << " in dimension " << Topology::dimension << std::endl;
  bool ret = true;

  LagrangeBasis<Topology,StorageField> basis(order);

  typedef Dune::LagrangePoints< Topology, StorageField > LagrangePoints;
  LagrangePoints points( order );

  std::vector< Dune::FieldVector< double, 1 > > y( basis.size() );

  for( unsigned int index = 0; index < points.size(); ++index )
  {
    basis.evaluate( points[ index ].point(), y );
    for( unsigned int i = 0; i < y.size(); ++i )
    {
      if( fabs( y[ i ] - double( i == index ) ) > 1e-10 )
      {
        std::cout << "ERROR: "
                  << index << " -> "
                  << "x = " << points[ index ].point()
                  << " (codim = " << points[ index ].localKey().codim() << ", "
                  << "subentity = " << points[ index ].localKey().subentity() << ", "
                  << "index = " << points[ index ].localKey().index() << "):"
                  << "    y[ " << i << " ] = " << y[ i ] << std::endl;
        ret = false;
      }
    }
  }
#if 0
  OrthonormalBasis< Topology, StorageField > basis( order );

  const unsigned int size = basis.size( order );
  std::vector< FieldVector< double, 1 > > y( size );

  std::vector< FieldVector< double, 1 > > m( size * size );
  for( unsigned int i = 0; i < size * size; ++i )
    m[ i ] = 0;

  GenericQuadrature< Topology > quadrature( 2*order+1 );
  const unsigned int quadratureSize = quadrature.size();
  for( unsigned int qi = 0; qi < quadratureSize; ++qi )
  {
    basis.evaluate( order, quadrature.point( qi ), y );
    for( unsigned int i = 0; i < size; ++i )
    {
      for( unsigned int j = 0; j < size; ++j )
        m[ i*size + j ] += quadrature.weight( qi ) * y[ i ] * y[ j ];
    }
  }

  for( unsigned int i = 0; i < size; ++i )
  {
    for( unsigned int j = 0; j < size; ++j )
    {
      const double value = m[ i*size + j ];
      if( fabs( value - double( i == j ) ) > 1e-10 ) {
        std::cout << "i = " << i << ", j = " << j << ": " << value << std::endl;
        ret = false;
      }
    }
  }
#endif
  if (!ret) {
    std::cout << "   FAILED !" << std::endl;
  }
  std::cout << std::endl;
  return ret;
}
int main ( int argc, char **argv )
{
  if( argc < 2 )
  {
    std::cerr << "Usage: " << argv[ 0 ] << " <p>" << std::endl;
    return 2;
  }

  const unsigned int order = atoi( argv[ 1 ] );
#ifdef TOPOLOGY
  return (test<TOPOLOGY>(order) ? 0 : 1 );
#else
  bool tests = true;
  tests &= test<Prism<Point> > (order);
  tests &= test<Pyramid<Point> > (order);

  tests &= test<Prism<Prism<Point> > > (order);
  tests &= test<Pyramid<Pyramid<Point> > >(order);

  tests &= test<Prism<Prism<Prism<Point> > > >(order);
  tests &= test<Prism<Pyramid<Pyramid<Point> > > >(order);
  tests &= test<Pyramid<Prism<Prism<Point> > > >(order);
  tests &= test<Pyramid<Pyramid<Pyramid<Point> > > >(order);

  tests &= test<Prism<Prism<Prism<Prism<Point> > > > >(order);
  tests &= test<Pyramid<Pyramid<Pyramid<Pyramid<Point> > > > >(order);
  return (tests ? 0 : 1);
#endif // TOPOLOGY
}
