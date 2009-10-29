// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <dune/finiteelements/generic/math/field.hh>
#include <dune/finiteelements/generic/orthonormalbasis/orthonormalbasis.hh>
#include <dune/finiteelements/generic/quadrature/genericquadrature.hh>

#if HAVE_ALGLIB
typedef amp::ampf< 128 > StorageField;
typedef amp::ampf< 512 > ComputeField;
#else
#if HAVE_GMP
typedef Dune::GMPField< 128 > StorageField;
typedef Dune::GMPField< 512 > ComputeField;
#else
typedef double StorageField;
typedef double ComputeField;
#endif
#endif

template <class Topology>
bool test(unsigned int order)
{
  bool ret = true;
  for (unsigned int o=order; o<=order; --o)
  {
    std::cout << "Testing " << Topology::name() << " in dimension " << Topology::dimension << " with order " << o << std::endl;
    typedef Dune::OrthonormalBasisFactory<Topology::dimension,StorageField,ComputeField> BasisFactory;
    const typename BasisFactory::Object &basis = *BasisFactory::template create<Topology>(o);

    const unsigned int size = basis.size( );

    std::vector< Dune::FieldVector< double, 1 > > y( size );

    std::vector< Dune::FieldVector< double, 1 > > m( size * size );
    for( unsigned int i = 0; i < size * size; ++i )
      m[ i ] = 0;

    Dune::GenericGeometry::GenericQuadrature< Topology, double > quadrature( 2*order+1 );
    const unsigned int quadratureSize = quadrature.size();
    for( unsigned int qi = 0; qi < quadratureSize; ++qi )
    {
      basis.evaluate( quadrature.point( qi ), y );
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
  }
  if (!ret) {
    std::cout << "   FAILED !" << std::endl;
  }
  std::cout << std::endl;
  return ret;
}
template <unsigned int dimension>
bool test(unsigned int topologyId, unsigned int order)
{
  bool ret = true;
  for (unsigned int o=order; o<=order; --o)
  {
    std::cout << "Testing " << topologyId << " in dimension " << dimension << " with order " << o << std::endl;
    typedef Dune::OrthonormalBasisFactory<dimension,StorageField,ComputeField> BasisFactory;
    const typename BasisFactory::Object &basis = *BasisFactory::create(topologyId,o);

    const unsigned int size = basis.size( );

    std::vector< Dune::FieldVector< double, 1 > > y( size );

    std::vector< Dune::FieldVector< double, 1 > > m( size * size );
    for( unsigned int i = 0; i < size * size; ++i )
      m[ i ] = 0;

    typedef typename Dune::GenericGeometry::GenericQuadratureProvider<dimension,double> QuadratureProvider;
    typedef typename QuadratureProvider::Quadrature Quadrature;
    const Quadrature &quadrature = QuadratureProvider::quadrature(topologyId,2*order+1);
    const unsigned int quadratureSize = quadrature.size();
    for( unsigned int qi = 0; qi < quadratureSize; ++qi )
    {
      basis.evaluate( quadrature.point( qi ), y );
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
  }
  if (!ret) {
    std::cout << "   FAILED !" << std::endl;
  }
  std::cout << std::endl;
  return ret;
}

int main ( int argc, char **argv )
{
  using namespace Dune;
  using namespace GenericGeometry;

  if( argc < 2 )
  {
    std::cerr << "Usage: " << argv[ 0 ] << " <p>" << std::endl;
    return 2;
  }

  const unsigned int order = atoi( argv[ 1 ] );
  bool tests = true;
#ifdef TOPOLOGY
  tests &= test<TOPOLOGY>(order);
  tests &= test<TOPOLOGY::dimension>(TOPOLOGY::id,order);
  return 0;
#else
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
