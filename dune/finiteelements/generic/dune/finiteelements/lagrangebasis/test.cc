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

  for( unsigned int index = 0; index < points.size(); ++index )
  {
    lagrange.evaluate( points[ index ].point(), y );
    std::cout << index << " -> "
              << "x = " << points[ index ].point()
              << " (codim = " << points[ index ].localKey().codim() << ", "
              << "subentity = " << points[ index ].localKey().subentity() << ", "
              << "index = " << points[ index ].localKey().index() << "):" << std::endl;

    for( unsigned int i = 0; i < y.size(); ++i )
    {
      if( fabs( y[ i ] - double( i == index ) ) > 1e-10 )
        std::cout << "    y[ " << i << " ] = " << y[ i ] << std::endl;
    }
  }
}
