// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <iostream>

#include <dune/alglib/multiprecision.hh>
#include <dune/alglib/matrix.hh>

#include <dune/finiteelements/lagrangebasis/lagrangepoints.hh>
#include <dune/finiteelements/lagrangebasis/interpolation.hh>
#include <dune/finiteelements/monomialbasis.hh>

#ifndef TOPOLOGY
#error "TOPOLOGY not defined."
#endif

using namespace Dune::GenericGeometry;

int main ( int argc, char **argv )
{
  if( argc < 2 )
  {
    std::cerr << "Usage: " << argv[ 0 ] << " <p>" << std::endl;
    return 1;
  }

  int p = atoi( argv[ 1 ] );

  typedef Dune::AlgLib::MultiPrecision< 256 > Field;
  typedef TOPOLOGY Topology;

  Dune::MonomialBasis< Topology, Field > basis( p );

  Dune::LocalLagrangeInterpolation< Topology, Field  > interpolation( p );

  Dune::AlgLib::Matrix< Field > matrix;
  interpolation.interpolate( basis, matrix );

  std::cout << "Matrix of evaluated base functions:" << std::endl;
  for( unsigned int row = 0; row < matrix.rows(); ++row )
  {
    for( unsigned int col = 0; col < matrix.cols(); ++col )
      std::cout << "   " << matrix( row, col );
    std::cout << std::endl;
  }
  std::cout << std::endl;

  matrix.invert();
  std::cout << "Inverse:" << std::endl;
  for( unsigned int row = 0; row < matrix.rows(); ++row )
  {
    for( unsigned int col = 0; col < matrix.cols(); ++col )
      std::cout << "   " << matrix( row, col );
    std::cout << std::endl;
  }
  std::cout << std::endl;
}
