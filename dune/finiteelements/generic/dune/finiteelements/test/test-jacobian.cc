// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <iostream>

#include <dune/finiteelements/lagrangepoints.hh>
#include <dune/finiteelements/monomialbasis.hh>
#include <dune/finiteelements/multiindex.hh>
#include <dune/alglib/multiprecision.hh>

#include <dune/finiteelements/lagrangebasis.hh>
#include <dune/finiteelements/quadrature/genericquadrature.hh>

#ifndef TOPOLOGY
#error "TOPOLOGY not defined."
#endif

using namespace Dune::GenericGeometry;
using namespace Dune;

int main ( int argc, char **argv )
{
  if( argc < 2 )
  {
    std::cerr << "Usage: " << argv[ 0 ] << " <p>" << std::endl;
    return 1;
  }

  int p = atoi( argv[ 1 ] );

  typedef TOPOLOGY Topology;
  const int dimension = Topology::dimension;
  typedef double StorageField;
  typedef amp::ampf< 256 > ComputeField;
  typedef double Field;

#if 0
  const int deriv = 1;
  typedef MultiIndex< dimension > MIF;
  typedef MonomialBasis< Topology, MIF > Basis;
  Basis basis;
  typedef MultiIndexEvaluator<Basis,Field> Evaluator;
#else
  const int deriv = 0;
  typedef LagrangeBasisProvider<dimension,StorageField,ComputeField> BasisProvider;
  typedef BasisProvider::Basis Basis;
  const Basis &basis = BasisProvider::basis(Topology::id,p);
  typedef StandardEvaluator<Basis> Evaluator;
#endif
  Evaluator eval(basis,p);
  typedef Evaluator::Iterator<deriv>::All Iterator;

  const unsigned int size = basis.size( );
  std::vector< Field > y( size );

  typedef LagrangePoints< Topology, Field > LagrangePoints;
  LagrangePoints points( p );

  std::cout << "Number of base functions:  " << size << std::endl;
  std::cout << "Number of Lagrange points: " << points.size() << std::endl;

  const LagrangePoints::iterator end = points.end();
  for( LagrangePoints::iterator it = points.begin(); it != end; ++it )
  {
    Iterator iter = eval.evaluateAll<deriv>(it->point());
    std::cout << "x = " << field_cast<double>(it->point())
              << " (codim = " << it->localKey().codim() << ", "
              << "subentity = " << it->localKey().subEntity() << ", "
              << "index = " << it->localKey().index() << "):" << std::endl;
    for (int i=0; !iter.done(); ++iter,++i)
    {
      // if (deriv==1)
      //   std::cout << "    y[ " << i << " ] = " << (*iter).jacobian;
      std::cout << "    y[ " << i << " ] = " << iter->value << std::endl;
    }
  }
}
