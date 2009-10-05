// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <iostream>
#include <vector>

#include <dune/alglib/multiprecision.hh>
#include <dune/finiteelements/tensor.hh>
const Dune::DerivativeLayout basisLayout = Dune::value;
const Dune::DerivativeLayout solutionLayout = Dune::value;

#include <dune/finiteelements/monomialbasis.hh>
#include <dune/finiteelements/multiindex.hh>
#include <dune/finiteelements/polynomialbasis.hh>
#include <dune/finiteelements/basisprint.hh>
#include <dune/finiteelements/basisevaluator.hh>

#ifndef TOPOLOGY
#error "TOPOLOGY not defined."
#endif

template <int dimR,int dimBasis>
struct TestMatrix;
template <int dimR>
struct TestMatrix<dimR,dimR>
{
  template <class Basis>
  TestMatrix(const Basis &basis)
    : size_(basis.size()) {}
  unsigned int colSize(int row) const {
    return size_*dimR;
  }
  unsigned int rowSize() const {
    return size_*dimR;
  }
  const double operator() ( int r, int c ) const
  {
    return (r+1)*(c+1);
    return pow(-1,c/2+r)*double(r)*double(c)*double(colSize(r)-c-1);
  }
  unsigned int size_;
};
template <int dimR>
struct TestMatrix<dimR,1>
{
  template <class Basis>
  TestMatrix(const Basis &basis)
    : size_(basis.size()) {}
  unsigned int colSize(int row) const {
    return size_;
  }
  unsigned int rowSize() const {
    return size_*dimR*dimR;
  }
  const double operator() ( int r, int c ) const
  {
    return (r/(dimR)+1)*(c*dimR+r%dimR+1);
    return pow(-1,c/2+r)*double(r)*double(c)*double(colSize(r)-c-1);
  }
  unsigned int size_;
};

using namespace Dune;
using namespace GenericGeometry;

template <class Topology>
void vecTest(unsigned int p)
{
  const int dimR = 2;
  const int dimBasis = 1;
  const int dimension = Topology::dimension;
  typedef MultiIndex< dimension > Field;

  typedef MonomialBasis< Topology, Field > Basis;
  Basis mbasis(p);
  typedef VectorialEvaluator<Basis,dimBasis,basisLayout> Evaluator;
  PolynomialBasisWithMatrix<Evaluator,SparseCoeffMatrix<double> > basis(mbasis);
  TestMatrix<dimR,dimBasis> matrix(mbasis);
  basis.fill(matrix);

  unsigned int size = basis.size()/2;
  std::cout << "Number of base functions:  " << size << std::endl;

  std::cout << ">>> Polynomial representation of the basis functions:" << std::endl;
  FieldVector< Field, dimension > x;
  for( int i = 0; i < dimension; ++i )
    x[ i ].set( i, 1 );

  std::cout << "Values: " << std::endl;
  std::vector< Derivatives<Field,dimension,dimR,0,solutionLayout> > val( size );
  for( unsigned int i = 0; i < val.size(); ++i )
    val[ i ] = -42.3456789;
  basis.template evaluate( x, val );
  std::cout << val << std::endl;

  std::cout << "Values+Jacobian: " << std::endl;
  std::vector< Derivatives<Field,dimension,dimR,1,solutionLayout> > deriv( size );
  for( unsigned int i = 0; i < deriv.size(); ++i )
    deriv[ i ] = -42.3456789;
  basis.template evaluate<1>( x, deriv );
  std::cout << deriv << std::endl;

  std::cout << "Values+Jacobian+Hessian: " << std::endl;
  std::vector< Derivatives<Field,dimension,dimR,2,solutionLayout> > hess( size );
  for( unsigned int i = 0; i < hess.size(); ++i )
    hess[ i ] = -42.3456789;
  basis.template evaluate<2>( x, hess );
  std::cout << hess  << std::endl;
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

  vecTest<Topology>(p);
}
