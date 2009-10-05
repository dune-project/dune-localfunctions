// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <dune/finiteelements/lagrangebasis.hh>
#include <dune/finiteelements/lagrangebasis/lagrangepoints.hh>
#include <dune/finiteelements/lagrangebasis/labattopoints.hh>
#include <dune/finiteelements/quadrature/genericquadrature.hh>
#include <dune/finiteelements/p13d/p13dlocalbasis.hh>
#include <dune/finiteelements/p23d/p23dlocalbasis.hh>
#include <dune/finiteelements/pk3d/pk3dlocalbasis.hh>

#define USE_GENERIC 1
const unsigned int iterations = 1;

using namespace Dune;
using namespace GenericGeometry;

struct SpecialBasis : public Pk3DLocalBasis<double,double,1>
                      // struct SpecialBasis :public P23DLocalBasis<double,double>
                      // struct SpecialBasis :public P13DLocalBasis<double,double>
{
  void evaluate(const Dune::FieldVector<double,3> &x,
                std::vector<Dune::FieldVector<double,1> > &ret)
  {
    evaluateFunction(x,ret);
  }

#if 0
  void evaluate(const Dune::FieldVector<double,2> &x,
                std::vector<Dune::FieldVector<double,1> > ret)
  {
    const double &b1 = x[0];
    const double &b2 = x[1];
    const double b0 = 1-b1-b2;
    ret[0] = b0*(2.*b0-1.) ;
    ret[1] = b1*(2.*b1-1.) ;
    ret[2] = b2*(2.*b2-1.) ;
    ret[3] = 4.*b1*b2;
    ret[4] = 4.*b2*b0;
    ret[5] = 4.*b0*b1;
    /*
       double bary[3]={1.-x[0]-x[1],x[0],x[1]};
       ret[0] = bary[0]*(2.*bary[0]-1.) ;
       ret[1] = bary[1]*(2.*bary[1]-1.) ;
       ret[2] = bary[2]*(2.*bary[2]-1.) ;
       ret[3] = 4.*bary[1]*bary[2];
       ret[4] = 4.*bary[2]*bary[0];
       ret[5] = 4.*bary[0]*bary[1];
     */
  }
#endif
};

template <class Topology>
bool test(unsigned int order, bool verbose = false) {

  typedef AlgLib::MultiPrecision<128> StorageField;
  // typedef double StorageField;
  typedef AlgLib::MultiPrecision<256> ComputeField;
  // typedef double ComputeField;

  bool ret = true;

  for (unsigned int o=1; o<=order; ++o)
  {
    std::cout << "# Testing " << Topology::name() << " in dimension " << Topology::dimension << " with order " << o << std::endl;

    typedef Dune::LagrangePoints< double, Topology::dimension > LagrangePoints;

    // typedef Dune::LagrangePointsCreator< StorageField, Topology::dimension > LagrangePointsCreator;
    typedef Dune::LobattoPointsCreator< double, Topology::dimension > LagrangePointsCreator;
    const LagrangePoints &points = LagrangePointsCreator::template lagrangePoints< Topology >( o );

#if USE_GENERIC
    // typedef LagrangeBasisProvider<Topology::dimension,StorageField,ComputeField> BasisProvider;
    typedef LobattoBasisProvider<Topology::dimension,StorageField,ComputeField> BasisProvider;
    const typename BasisProvider::Basis &basis = BasisProvider::basis(Topology::id,o);
    std::vector< Dune::FieldVector< double, 1 > > y( basis.size() );
    for (unsigned int count = 0; count < iterations; ++count)
    {
      for( unsigned int index = 0; index < points.size(); ++index )
      {
        if (verbose)
          std::cout << index << "   " << points[ index ].point() << " "
                    << points[ index ].localKey()
                    << std::endl;
        basis.evaluate( points[ index ].point(), y );
        bool first = true;
        if (iterations==1)
        {
          for( unsigned int i = 0; i < y.size(); ++i )
          {
            if( fabs( y[ i ] - double( i == index ) ) > 1e-10 )
            {
              if (first && verbose) {
                std::cout << "ERROR: "
                          << index << " -> "
                          << "x = " << points[ index ].point()
                          << " (codim = " << points[ index ].localKey().codim() << ", "
                          << "subentity = " << points[ index ].localKey().subEntity() << ", "
                          << "index = " << points[ index ].localKey().index() << "):" << std::endl;
                first = false;
              }
              if (verbose)
                std::cout << "         y[ " << i << " ] = " << y[ i ] << " "
                          << "         error : " << fabs( y[ i ] - double( i == index ) )
                          << std::endl;
              ret = false;
            }
          }
        }
      }
    }
#else
    SpecialBasis specialBasis;
    std::vector< Dune::FieldVector< double, 1 > > y( points.size() );
    for (unsigned int count = 0; count < iterations; ++count)
    {
      for( unsigned int index = 0; index < points.size(); ++index )
      {
        specialBasis.evaluate( points[ index ].point(), y );
        bool first = true;
        if (iterations==1)
        {
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
              std::cout << "         y[ " << i << " ] = " << y[ i ] << std::endl;
              ret = false;
            }
          }
        }
      }
    }
#endif
    if (verbose)
      std::cout << std::endl << std::endl << std::endl;
  }
  if (!ret) {
    std::cout << "   FAILED !" << std::endl;
  }
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
  return (test<TOPOLOGY>(order,true) ? 0 : 1 );
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
