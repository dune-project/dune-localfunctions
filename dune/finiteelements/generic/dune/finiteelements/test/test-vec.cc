// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <iostream>
#include <vector>

#include <dune/alglib/multiprecision.hh>
#include <dune/finiteelements/tensor.hh>

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

template <class Topology,
    int dimBasis,
    DerivativeLayout basisLayout,
    int dimR,
    DerivativeLayout solutionLayout>
void vecTest(int testNr,unsigned int p)
{
  std::cout << "Starting on test : " << testNr << std::endl;
  std::stringstream name;
  name << "vectest-" << testNr << ".out";
  std::ofstream out(name.str().c_str());
  const int dimension = Topology::dimension;
  typedef MultiIndex< dimension > Field;

  typedef MonomialBasis< Topology, Field > Basis;
  Basis mbasis(p);
  typedef VectorialEvaluator<Basis,dimBasis,basisLayout> Evaluator;
  const unsigned int blockSize = (dimBasis==dimR) ? 1 : dimR;
  PolynomialBasisWithMatrix<Evaluator,SparseCoeffMatrix<double,blockSize> > basis(mbasis);
  TestMatrix<dimR,dimBasis> matrix(mbasis);
  basis.fill(matrix);

  unsigned int size = basis.size();
  out << "Number of base functions:  " << size << std::endl;

  out << ">>> Polynomial representation of the basis functions:" << std::endl;
  FieldVector< Field, dimension > x;
  for( int i = 0; i < dimension; ++i )
    x[ i ].set( i, 1 );

  out << "Values: " << std::endl;
  std::vector< Derivatives<Field,dimension,dimR,0,solutionLayout> > val( size );
  for( unsigned int i = 0; i < val.size(); ++i )
    val[ i ] = -42.3456789;
  basis.template evaluate( x, val );
  out << val << std::endl;

  out << "Values+Jacobian: " << std::endl;
  std::vector< Derivatives<Field,dimension,dimR,1,solutionLayout> > deriv( size );
  for( unsigned int i = 0; i < deriv.size(); ++i )
    deriv[ i ] = -42.3456789;
  basis.template evaluate<1>( x, deriv );
  out << deriv << std::endl;

  out << "Values+Jacobian+Hessian: " << std::endl;
  std::vector< Derivatives<Field,dimension,dimR,2,solutionLayout> > hess( size );
  for( unsigned int i = 0; i < hess.size(); ++i )
    hess[ i ] = -42.3456789;
  basis.template evaluate<2>( x, hess );
  out << hess  << std::endl;
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

  // problem: vecTest<Topology,1,value,2,value>(1,p);
  vecTest<Topology,1,derivative,2,derivative>(2,p);
  vecTest<Topology,2,value,2,value>(3,p);
  vecTest<Topology,2,derivative,2,derivative>(4,p);
  // problem: vecTest<Topology,1,derivative,2,value>(5,p);
  vecTest<Topology,1,value,2,derivative>(6,p);
  // missing: vecTest<Topology,2,value,2,derivative>(7,p);
  vecTest<Topology,2,derivative,2,value>(8,p);
}
