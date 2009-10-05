// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <dune/finiteelements/lagrangebasis.hh>
#include <dune/finiteelements/quadrature/genericquadrature.hh>

using namespace Dune;
using namespace GenericGeometry;

template <class Topology>
bool test(unsigned int order) {
  typedef AlgLib::MultiPrecision<128> StorageField;
  // typedef double StorageField;
  typedef AlgLib::MultiPrecision<512> ComputeField;
  // typedef double ComputeField;

  bool ret = true;

  for (unsigned int o=1; o<=order; ++o)
  {
    std::cout << "Testing " << Topology::name() << " in dimension " << Topology::dimension << " with order " << o << std::endl;
    typedef LagrangeBasisProvider<Topology::dimension,StorageField> BasisProvider;
    const typename BasisProvider::Basis &basis = BasisProvider::basis(Topology::id,o);

    typedef Dune::LagrangePoints< Topology, StorageField > LagrangePoints;
    LagrangePoints points( o );

    std::vector< Dune::FieldVector< double, 1 > > y( basis.size() );
    for( unsigned int index = 0; index < points.size(); ++index )
    {
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
                      << "subentity = " << points[ index ].localKey().subentity() << ", "
                      << "index = " << points[ index ].localKey().index() << "):" << std::endl;
            first = false;
          }
          std::cout << "         y[ " << i << " ] = " << y[ i ] << std::endl;
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
  if( argc < 2 )
  {
    std::cerr << "Usage: " << argv[ 0 ] << " <p>" << std::endl;
    return 2;
  }

  const unsigned int order = atoi( argv[ 1 ] );
#ifdef TOPOLOGY
  return (test<TOPOLOGY>(order,fixed) ? 0 : 1 );
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
