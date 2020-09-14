// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <dune/localfunctions/utility/field.hh>
#include <dune/localfunctions/utility/basisprint.hh>

#include <dune/localfunctions/lagrange/equidistantpoints.hh>
#include <dune/localfunctions/lagrange/lagrangebasis.hh>

/**
 * \file
 * \brief Performs some tests for the generic Lagrange
 *        shape functions on simplices.
 *
 * The topology can be chosen at compile time by setting TOPOLOGY
 * to a string like
 * \code
 * Pyramid<Pyramid<Point> > >
 * \endcode
 * which generates a 2d simplex. If TOPOLOGY is not set, all
 * topologies up to 4d are tested. Note, this may lead to prolonged
 * compiler runs.
 *
 * For debugging purpuse the functions and the derivatives can be
 * printed. You have to define the macro TEST_OUTPUT_FUNCTIONS to
 * activate this function.
 */

#if HAVE_GMP
typedef Dune::GMPField< 128 > StorageField;
typedef Dune::GMPField< 512 > ComputeField;
#else
typedef double StorageField;
typedef double ComputeField;
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

  bool ret = true;

  for (unsigned int o = 0; o <= order; ++o)
  {
    const typename LagrangeCoefficientsFactory::Object *pointsPtr = LagrangeCoefficientsFactory::template create< Topology >( o );

    if ( pointsPtr == 0)
      continue;

    std::cout << "# Testing " << Topology::name() << " in dimension " << Topology::dimension << " with order " << o << std::endl;

    typename BasisFactory::Object &basis = *BasisFactory::template create<Topology>(o);

    ret |= test(basis,*pointsPtr,verbose);

    // define the macro TEST_OUTPUT_FUNCTIONS to output files containing functions and
    // derivatives in a human readabible form (aka LaTeX source)
#ifdef TEST_OUTPUT_FUNCTIONS
    std::stringstream name;
    name << "lagrange_" << Topology::name() << "_p" << o << ".basis";
    std::ofstream out(name.str().c_str());
    Dune::basisPrint<0,BasisFactory,typename BasisFactory::StorageField,Topology>(out,basis);
    Dune::basisPrint<1,BasisFactory,typename BasisFactory::StorageField,Topology>(out,basis);
#endif // TEST_OUTPUT_FUNCTIONS

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

#ifdef CHECKDIM
  #if CHECKDIM==1
      #define CHECKDIM1
  #elif CHECKDIM==2
      #define CHECKDIM2
  #elif CHECKDIM==3
      #define CHECKDIM3
  #elif CHECKDIM==4
      #define CHECKDIM4
  #endif
#else
  #define CHECKDIM1
  #define CHECKDIM2
  #define CHECKDIM3
  #define CHECKDIM4
#endif



int main ( int argc, char **argv )
{
  using namespace Dune;
  using namespace Impl;

  unsigned int order = (argc < 2) ? 5 : atoi(argv[1]);

  if (argc < 2)
  {
    std::cerr << "Usage: " << argv[ 0 ] << " <p>" << std::endl
              << "Using default order of " << order << std::endl;
  }
#ifdef TOPOLOGY
  return (test<TOPOLOGY>(order) ? 0 : 1 );
#else
  bool tests = true;

#ifdef CHECKDIM1
  tests &= test<Prism<Point> > (order);
  tests &= test<Pyramid<Point> > (order);
#endif

#ifdef CHECKDIM2
  tests &= test<Prism<Prism<Point> > > (order);
  tests &= test<Pyramid<Pyramid<Point> > >(order);
#endif

#ifdef CHECKDIM3
  tests &= test<Prism<Prism<Prism<Point> > > >(order);
  tests &= test<Prism<Pyramid<Pyramid<Point> > > >(order);
  tests &= test<Pyramid<Pyramid<Pyramid<Point> > > >(order);
#endif

  // tests &= test<Pyramid<Prism<Prism<Point> > > >(order);
  std::cout << "NOT CHECKING PYRAMID!" << std::endl;

  // reduce tested order to 4 in 4d unless explicitly asked for more
  if (argc < 2)
    order = 4;

#ifdef CHECKDIM4
  tests &= test<Prism<Prism<Prism<Prism<Point> > > > >(order);
  tests &= test<Pyramid<Pyramid<Pyramid<Pyramid<Point> > > > >(order);
#endif

  return (tests ? 0 : 1);
#endif // TOPOLOGY
}
