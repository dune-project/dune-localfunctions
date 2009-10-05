// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <dune/finiteelements/raviartthomas/raviartthomasbasis.hh>

using namespace Dune;
using namespace GenericGeometry;

template <class Topology>
bool test(unsigned int order) {
  // typedef AlgLib::MultiPrecision<128> StorageField;
  typedef double StorageField;
  typedef AlgLib::MultiPrecision<512> ComputeField;
  // typedef double ComputeField;

  bool ret = true;

  for (unsigned int o=order; o<=order; ++o)
  {
    std::cout << "Testing " << Topology::name() << " in dimension " << Topology::dimension << " with order " << o << std::endl;
    typedef RaviartThomasBasisProvider<Topology::dimension,StorageField,ComputeField> BasisProvider;
    const typename BasisProvider::Basis &basis = BasisProvider::basis(Topology::id,o);

    typedef Dune::LagrangePoints< Topology, StorageField > LagrangePoints;
    LagrangePoints points( 1 );

    std::vector< Dune::FieldVector< double, Topology::dimension > > y( basis.size() );
    for( unsigned int index = 0; index < points.size(); ++index )
    {
      basis.evaluate( points[ index ].point(), y );
      bool first = true;
      std::cout << "At point points[ " << index << " ] = " << points[ index ].point() << std::endl;
      for( unsigned int i = 0; i < y.size(); ++i)
      {
        std::cout << "  y [ " << i << " ] = " << y[i] << std::endl;
        first = false;
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
  return (test<TOPOLOGY>(order) ? 0 : 1 );
#else
  bool tests = true;
  tests &= test<Pyramid<Point> > (order);

  tests &= test<Pyramid<Pyramid<Point> > >(order);

  tests &= test<Pyramid<Pyramid<Pyramid<Point> > > >(order);

  tests &= test<Pyramid<Pyramid<Pyramid<Pyramid<Point> > > > >(order);
  return (tests ? 0 : 1);
#endif // TOPOLOGY
}
