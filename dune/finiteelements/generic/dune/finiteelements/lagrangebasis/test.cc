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

template <class Basis,class Points>
bool test(const Basis &basis, const Points &points, bool verbose)
{
  bool ret = true;
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
  return ret;
}

template <class Topology>
bool test(unsigned int order, bool verbose = false)
{
  typedef Dune::LagrangeBasisFactory<Dune::EquidistantPointSet,Topology::dimension,StorageField,ComputeField> BasisFactory;
  typedef Dune::LagrangeCoefficientsFactory< Dune::EquidistantPointSet,  Topology::dimension,double > LagrangeCoefficientsFactory;
  // typedef Dune::LagrangeBasisFactory<Dune::LobattoPointSet,Topology::dimension,StorageField,ComputeField> BasisFactory;
  // typedef Dune::LagrangeCoefficientsFactory< Dune::LobattoPointSet, Topology::dimension, double > LagrangeCoefficientsFactory;

  bool ret = true;

  for (unsigned int o=1; o<=order; ++o)
  {
    const typename LagrangeCoefficientsFactory::Object *pointsPtr = LagrangeCoefficientsFactory::template create< Topology >( o );

    if ( pointsPtr == 0)
      continue;

    std::cout << "# Testing " << Topology::name() << " in dimension " << Topology::dimension << " with order " << o << std::endl;

    typename BasisFactory::Object &basis = *BasisFactory::template create<Topology>(o);

    ret |= test(basis,*pointsPtr,verbose);

    LagrangeCoefficientsFactory::release( pointsPtr );
    BasisFactory::release( &basis );
  }

  if (verbose)
    std::cout << std::endl << std::endl << std::endl;
  if (!ret) {
    std::cout << "   FAILED !" << std::endl;
  }
  return ret;
}
template <unsigned int dimension>
bool test(unsigned int topologyId, unsigned int order, bool verbose = false)
{
  typedef Dune::LagrangeBasisFactory<Dune::EquidistantPointSet,dimension,StorageField,ComputeField> BasisFactory;
  typedef Dune::LagrangeCoefficientsFactory< Dune::EquidistantPointSet,  dimension,double > LagrangeCoefficientsFactory;
  // typedef Dune::LagrangeBasisFactory<Dune::LobattoPointSet,Topology::dimension,StorageField,ComputeField> BasisFactory;
  // typedef Dune::LagrangeCoefficientsFactory< Dune::LobattoPointSet, Topology::dimension, double > LagrangeCoefficientsFactory;

  bool ret = true;

  for (unsigned int o=1; o<=order; ++o)
  {
    const typename LagrangeCoefficientsFactory::Object *pointsPtr = LagrangeCoefficientsFactory::create( topologyId, o );

    if ( pointsPtr == 0)
      continue;

    std::cout << "# Testing " << topologyId << " in dimension " << dimension << " with order " << o << std::endl;

    typename BasisFactory::Object &basis = *BasisFactory::create( topologyId, o );

    ret |= test(basis,*pointsPtr,verbose);

    LagrangeCoefficientsFactory::release( pointsPtr );
    BasisFactory::release( &basis );
  }

  if (verbose)
    std::cout << std::endl << std::endl << std::endl;
  if (!ret) {
    std::cout << "   FAILED !" << std::endl;
  }
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
