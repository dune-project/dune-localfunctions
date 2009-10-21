// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <dune/finiteelements/common/field.hh>

#include <dune/finiteelements/lagrangebasis/equidistantpoints.hh>
#if HAVE_ALGLIB
#include <dune/finiteelements/lagrangebasis/lobattopoints.hh>
#endif
#include <dune/finiteelements/lagrangebasis/lagrangebasis.hh>
#include <dune/finiteelements/quadrature/genericquadrature.hh>

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
bool test(unsigned int order, bool verbose = false) {

  bool ret = true;

  for (unsigned int o=1; o<=order; ++o)
  {
    std::cout << "# Testing " << Topology::name() << " in dimension " << Topology::dimension << " with order " << o << std::endl;

    // typedef Dune::LagrangeCoefficientsFactory< Dune::EquidistantCoefficients,  Topology::dimension,double > LagrangeCoefficientsFactory;
    typedef Dune::LagrangeCoefficientsFactory< Dune::LobattoCoefficients, Topology::dimension, double > LagrangeCoefficientsFactory;
    const typename LagrangeCoefficientsFactory::Object *pointsPtr = LagrangeCoefficientsFactory::template create< Topology >( o );

    if ( pointsPtr == 0)
      continue;
    const typename LagrangeCoefficientsFactory::Object &points = *pointsPtr;

    // typedef Dune::LagrangeBasisFactory<Dune::EquidistantCoefficients,Topology::dimension,StorageField,ComputeField> BasisFactory;
    typedef Dune::LagrangeBasisFactory<Dune::LobattoCoefficients,Topology::dimension,StorageField,ComputeField> BasisFactory;
    typename BasisFactory::Object &basis = *BasisFactory::template create<Topology>(o);

    std::vector< Dune::FieldVector< double, 1 > > y( basis.size() );
    for( unsigned int index = 0; index < points.size(); ++index )
    {
      if (verbose)
        std::cout << index << "   " << points[ index ].point() << " "
                  << points[ index ].localKey()
                  << std::endl;
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
                      << "subentity = " << points[ index ].localKey().subEntity() << ", "
                      << "index = " << points[ index ].localKey().index() << "):" << std::endl;
            first = false;
          }
          if (1)
            std::cout << "         y[ " << i << " ] = " << y[ i ] << " "
                      << "         error : " << fabs( y[ i ] - double( i == index ) )
                      << std::endl;
          ret = false;
        }
      }
    }

    // add release(basis), release(points)

    if (verbose)
      std::cout << std::endl << std::endl << std::endl;
  }
  if (!ret) {
    std::cout << "   FAILED !" << std::endl;
  }
  return ret;
}
template <unsigned int dimension>
bool test(unsigned int topologyId, unsigned int order, bool verbose = false)
{
#if 0
  bool ret = true;

  for (unsigned int o=1; o<=order; ++o)
  {
    std::cout << "# Testing " << topologyId << " in dimension " << dimension << " with order " << o << std::endl;

    // typedef Dune::LagrangeCoefficientsFactory< double, dimension > LagrangeCoefficientsFactory;
    typedef Dune::LobattoCoefficientsFactory< double, dimension > LagrangeCoefficientsFactory;
    const typename LagrangeCoefficientsFactory::LagrangeCoefficients *pointsPtr = LagrangeCoefficientsFactory::create( topologyId, o );

    if ( pointsPtr == 0)
      continue;
    const typename LagrangeCoefficientsFactory::LagrangeCoefficients &points = *pointsPtr;

    // typedef Dune::LagrangeBasisFactory<dimension,Dune::LagrangeCoefficientsFactory,StorageField,ComputeField> BasisFactory;
    typedef Dune::LagrangeBasisFactory<dimension,Dune::LobattoCoefficientsFactory,StorageField,ComputeField> BasisFactory;
    typename BasisFactory::Object &basis = *BasisFactory::create(topologyId,o);

    std::vector< Dune::FieldVector< double, 1 > > y( basis.size() );
    for( unsigned int index = 0; index < points.size(); ++index )
    {
      if (verbose)
        std::cout << index << "   " << points[ index ].point() << " "
                  << points[ index ].localKey()
                  << std::endl;
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
                      << "subentity = " << points[ index ].localKey().subEntity() << ", "
                      << "index = " << points[ index ].localKey().index() << "):" << std::endl;
            first = false;
          }
          if (1)
            std::cout << "         y[ " << i << " ] = " << y[ i ] << " "
                      << "         error : " << fabs( y[ i ] - double( i == index ) )
                      << std::endl;
          ret = false;
        }
      }
    }

    // add release(basis), release(points)

    if (verbose)
      std::cout << std::endl << std::endl << std::endl;
  }
  if (!ret) {
    std::cout << "   FAILED !" << std::endl;
  }
  return ret;
#endif
  return 1;
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
#ifdef TOPOLOGY
  test<TOPOLOGY>(order,false);
  test<TOPOLOGY::dimension>(TOPOLOGY::id,order,false);
#else
  bool tests = true;
  tests &= test<Prism<Point> > (order);
  tests &= test<Pyramid<Point> > (order);

  tests &= test<Prism<Prism<Point> > > (order);
  tests &= test<Pyramid<Pyramid<Point> > >(order);

  tests &= test<Prism<Prism<Prism<Point> > > >(order);
  tests &= test<Prism<Pyramid<Pyramid<Point> > > >(order);
  tests &= test<Pyramid<Pyramid<Pyramid<Point> > > >(order);

  // tests &= test<Pyramid<Prism<Prism<Point> > > >(order);
  std::cout << "NOT CHECKING PYRAMID!" << std::endl;

  tests &= test<Prism<Prism<Prism<Prism<Point> > > > >(order);
  tests &= test<Pyramid<Pyramid<Pyramid<Pyramid<Point> > > > >(order);
  return (tests ? 0 : 1);
#endif // TOPOLOGY
}
