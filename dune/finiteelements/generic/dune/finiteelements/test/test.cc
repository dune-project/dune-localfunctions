// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <iostream>
#include <vector>

#include <dune/alglib/multiprecision.hh>
#include <dune/finiteelements/quadrature/genericquadrature.hh>
#include <dune/finiteelements/lagrangepoints.hh>
#include <dune/finiteelements/monomialbasis.hh>
#include <dune/finiteelements/multiindex.hh>
#include <dune/finiteelements/basisevaluator.hh>
#include <dune/finiteelements/polynomialbasis.hh>

#ifndef TOPOLOGY
#error "TOPOLOGY not defined."
#endif

template <int dimR>
struct TestMatrix
{
  template <class Basis>
  TestMatrix(const Basis &basis)
    : size_(basis.size()) {}
  int colSize(int row) const {
    return size_*dimR;
  }
  int rowSize() const {
    return size_*dimR;
  }
  const double operator() ( int r, int c ) const
  {
    return double(r+1)*double(c+1);
  }
  unsigned int size_;
};

using namespace Dune;
using namespace GenericGeometry;

template <class Topology>
void lagrangePointTest(unsigned int p)
{
  typedef double Field;

  typedef MonomialBasis< Topology, Field > Basis;
  Basis basis(p);

  unsigned int size = basis.size();
  std::cout << "Number of base functions:  " << size << std::endl;

  typedef LagrangePoints< Topology, Field > LagrangePoints;
  LagrangePoints points( p );
  std::cout << "Number of Lagrange points: " << points.size() << std::endl;
  std::vector< Field > y( size );

  std::cout << "Testing monomial basis: " << std::endl;
  const typename LagrangePoints::iterator end = points.end();
  for( typename LagrangePoints::iterator it = points.begin(); it != end; ++it )
  {
    basis.evaluate( it->point(), &(y[0]) );
    std::cout << "x = " << field_cast<double>(it->point())
              << " (codim = " << it->localKey().codim() << ", "
              << "subentity = " << it->localKey().subEntity() << ", "
              << "index = " << it->localKey().index() << "):" << std::endl;
    for( unsigned int i = 0; i < size; ++i )
      std::cout << "    y[ " << i << " ] = " << field_cast<double>(y[ i ]) << std::endl;
  }

}
template <class Topology>
void quadratureTest(unsigned int p)
{
  typedef double Field;

  typedef MonomialBasis< Topology, Field > Basis;
  Basis basis(p);

  unsigned int size = basis.size();
  std::cout << "Number of base functions:  " << size << std::endl;
  std::cout << std::endl;
  std::cout << ">>> Testing quadrature of order " << (2*p+1) << "..." << std::endl;

  std::vector< FieldVector< Field, 1 > > yquad( size );
  for( unsigned int i = 0; i < size; ++i )
    yquad[ i ] = 0;

  std::vector< Field > y( size );
  GenericQuadrature< Topology,Field > quadrature( 2*p+1 );
  const unsigned int quadratureSize = quadrature.size();
  for( unsigned int qi = 0; qi < quadratureSize; ++qi )
  {
    basis.evaluate( quadrature.point( qi ), y );
    std::cout << "x = " << quadrature.point( qi ) << ":" << std::endl;
    for( unsigned int i = 0; i < size; ++i )
    {
      yquad[ i ] += quadrature.weight( qi ) * y[ i ];
      std::cout << "    y[ " << i << " ] = " << y[ i ] << std::endl;
    }
  }

  std::vector< double > yint( size );
  basis.integral( yint );
  for( unsigned int i = 0; i < size; ++i )
  {
    if( fabs( yquad[ i ] - yint[ i ] ) < 1e-10 )
      continue;
    std::cerr << "Quadrature and Integral differ for basis function " << i << "." << std::endl;
    std::cout << "    quadrature: " << yquad[ i ] << std::endl;
    std::cout << "    integral:   " << yint[ i ] << std::endl;
  }
}
template <class Topology>
void polynomialBaseTest(unsigned int p)
{
  typedef double Field;

  typedef MonomialBasis< Topology, Field > Basis;
  Basis basis(p);

  unsigned int size = basis.size();
  std::cout << "Number of base functions:  " << size << std::endl;

  typedef LagrangePoints< Topology, Field > LagrangePoints;
  LagrangePoints points( p );
  std::cout << "Number of Lagrange points: " << points.size() << std::endl;
  std::vector< Field > y( size );

  std::cout << "Testing polynomial basis: " << std::endl;
  typedef StandardEvaluator<Basis> Evaluator;
  Evaluator eval(basis);
  PolynomialBasisWithMatrix<Evaluator> pBasis(basis);
  TestMatrix<1> matrix(basis);
  pBasis.fill(matrix);
  const typename LagrangePoints::iterator end = points.end();
  for( typename LagrangePoints::iterator it = points.begin(); it != end; ++it )
  {
    pBasis.evaluate( it->point(), y );
    std::cout << "x = " << field_cast<double>(it->point())
              << " (codim = " << it->localKey().codim() << ", "
              << "subentity = " << it->localKey().subEntity() << ", "
              << "index = " << it->localKey().index() << "):" << std::endl;
    for( unsigned int i = 0; i < size; ++i )
      std::cout << "    y[ " << i << " ] = " << field_cast<double>(y[ i ]) << std::endl;
  }
}
template <class Topology>
void multiIndexTest(unsigned int p)
{
  const int dimension = Topology::dimension;
  typedef MultiIndex< dimension > Field;

  typedef MonomialBasis< Topology, Field > Basis;
  Basis basis(p);

  unsigned int size = basis.size();
  std::cout << "Number of base functions:  " << size << std::endl;

  std::cout << ">>> Polynomial representation of the basis functions:" << std::endl;
  FieldVector< Field, dimension > x;
  for( int i = 0; i < dimension; ++i )
    x[ i ].set( i, 1 );
  std::vector< Field > val( size );
  basis.evaluate( x, val );
  std::cout << val << std::endl;

  std::cout << ">>> integral of basis functions:" << std::endl;
  std::vector< Field > integral( size );
  basis.integral( integral );
  std::cout << integral << std::endl;

  const int dimR = 2;
  std::cout << ">>> Polynomial basis: " << std::endl;
  typedef VectorialEvaluator<Basis,dimR> Evaluator;
  Evaluator eval(basis);
  PolynomialBasisWithMatrix<Evaluator,SparseCoeffMatrix<double> > pBasis(basis);
  TestMatrix<dimR> matrix(basis);
  pBasis.fill(matrix);
  std::vector< Dune::FieldVector<Field,dimR> > y( pBasis.size() );
  pBasis.evaluate( x, y );
  std::cout << y << std::endl;
  //for( unsigned int i = 0; i < size; ++i )
  //  std::cout << "    y[ " << i << " ] = " << y[ i ] << std::endl;
}

int main ( int argc, char **argv )
{
  typedef TOPOLOGY Topology;
  if( argc < 2 )
  {
    std::cerr << "Usage: " << argv[ 0 ] << " <p>" << std::endl;
    return 1;
  }
  int p = atoi( argv[ 1 ] );

  // lagrangePointTest<Topology>(p);
  // quadratureTest<Topology>(p);
  // polynomialBaseTest<Topology>(p);
  multiIndexTest<Topology>(p);
}
