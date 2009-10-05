// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include "lagrangebasis.hh"
using namespace Dune;
using namespace GenericGeometry;
int main ( int argc, char **argv )
{
  int order = atoi(argv[1]);
  LagrangeBasis<TOPOLOGY,double> lagrange(order);

  typedef Dune::LagrangePoints< TOPOLOGY, double > LagrangePoints;
  LagrangePoints points( order );

  std::vector< Dune::FieldVector< double, 1 > > y( lagrange.size() );

  int index = 0;
  const LagrangePoints::iterator end = points.end();
  for( LagrangePoints::iterator it = points.begin(); it != end; ++it,++index )
  {
    lagrange.evaluate( it->point(), y );
    std::cout << index << " -> "
              << "x = " << it->point()
              << " (codim = " << it->localKey().codim() << ", "
              << "subentity = " << it->localKey().subentity() << ", "
              << "index = " << it->localKey().index() << "):" << std::endl;
    for( unsigned int i = 0; i < y.size(); ++i ) {
      if ( i  == index ) {
        if (fabs(y[ i ] - 1.0) > 1e-10)
          std::cout << "    y[ " << i << " ] = " << y[ i ] << std::endl;
      } else {
        if (fabs(y[ i ]) > 1e-10)
          std::cout << "    y[ " << i << " ] = " << y[ i ] << std::endl;
      }
    }
  }
}
