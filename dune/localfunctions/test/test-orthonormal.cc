// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <dune/geometry/type.hh>
#include <dune/geometry/quadraturerules.hh>

#include <dune/localfunctions/utility/field.hh>
#include <dune/localfunctions/utility/basisprint.hh>
#include <dune/localfunctions/orthonormal/orthonormalbasis.hh>

/**
 * \file
 * \brief Performs some tests for the generic orthonormal
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

template <class Topology>
bool test(unsigned int order)
{
  bool ret = true;
  Dune::GeometryType gt(Topology::id, Topology::dimension);
  for (unsigned int o = 0; o <= order; ++o)
  {
    std::cout << "Testing " << Topology::name() << " in dimension " << Topology::dimension << " with order " << o << std::endl;
    typedef Dune::OrthonormalBasisFactory<Topology::dimension,StorageField,ComputeField> BasisFactory;
    const typename BasisFactory::Object &basis = *BasisFactory::template create<Topology>(o);

    const unsigned int size = basis.size( );

    std::vector< Dune::FieldVector< double, 1 > > y( size );

    std::vector< Dune::FieldVector< double, 1 > > m( size * size );
    for( unsigned int i = 0; i < size * size; ++i )
      m[ i ] = 0;

    const Dune::QuadratureRule<double,Topology::dimension> &quadrature =
      Dune::QuadratureRules<double,Topology::dimension>::rule(gt,2*order+1);
    const unsigned int quadratureSize = quadrature.size();
    for( unsigned int qi = 0; qi < quadratureSize; ++qi )
    {
      basis.evaluate( quadrature[qi].position(), y );
      for( unsigned int i = 0; i < size; ++i )
      {
        for( unsigned int j = 0; j < size; ++j )
          m[ i*size + j ] += quadrature[qi].weight() * y[ i ] * y[ j ];
      }
    }

    for( unsigned int i = 0; i < size; ++i )
    {
      for( unsigned int j = 0; j < size; ++j )
      {
        const double value = m[ i*size + j ];
        if( std::abs( value - double( i == j ) ) > 1200.*Dune::Zero<double>::epsilon() ) {
          std::cout << "i = " << i << ", j = " << j << ": " << std::abs( value - double( i == j ) ) << std::endl;
          ret = false;
        }
      }
    }

    // define the macro TEST_OUTPUT_FUNCTIONS to output files containing functions and
    // derivatives in a human readabible form (aka LaTeX source)
#ifdef TEST_OUTPUT_FUNCTIONS
    std::stringstream name;
    name << "orthonormal_" << Topology::name() << "_p" << o << ".basis";
    std::ofstream out(name.str().c_str());
    Dune::basisPrint<0,BasisFactory,typename BasisFactory::StorageField,Topology>(out,basis);
    Dune::basisPrint<1,BasisFactory,typename BasisFactory::StorageField,Topology>(out,basis);
#endif // TEST_OUTPUT_FUNCTIONS

    BasisFactory::release(&basis);
  }
  if (!ret) {
    std::cout << "   FAILED !" << std::endl;
  }
  std::cout << std::endl;
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

  const unsigned int order = (argc < 2) ? 5 : atoi(argv[1]);

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
  tests &= test<Pyramid<Prism<Prism<Point> > > >(order);
  tests &= test<Pyramid<Pyramid<Pyramid<Point> > > >(order);
#endif

#ifdef CHECKDIM4
  tests &= test<Prism<Prism<Prism<Prism<Point> > > > >(order);
  tests &= test<Pyramid<Pyramid<Pyramid<Pyramid<Point> > > > >(order);
#endif

  return (tests ? 0 : 1);
#endif // TOPOLOGY
}
