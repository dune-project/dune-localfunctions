// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include "test-space.hh"

#ifdef DGLAGRANGE
#include <dune/finiteelements/lagrangebasis/dgspace.hh>
typedef Dune::LagrangeDGSpace< GridView, double > Space;
std::string name("dglagrange");
#endif
#ifdef LAGRANGE
#include <dune/finiteelements/lagrangebasis/space.hh>
typedef Dune::LagrangeSpace< GridView, double > Space;
std::string name("lagrange");
#endif
#ifdef LOBATTO
#include <dune/finiteelements/lagrangebasis/lobattospace.hh>
typedef Dune::LobattoLagrangeSpace< GridView, double > Space;
std::string name("lobatto");
#endif
#ifdef ONB
#include <dune/finiteelements/orthonormalbasis/dgspace.hh>
typedef Dune::OrthonormalDGSpace< GridView, double > Space;
std::string name("onb");
#endif
#ifdef MONOMIALS
#endif

typedef Dune::Interpolation< Space > Interpolation;
typedef Dune::VTKFunctionWrapper< Space > VTKFunction;

int main ( int argc, char **argv )
{
  if( argc < 5 )
  {
    std::cerr << "Usage: " << argv[ 0 ] << " <dgf-file> <order> <level> <problem> " << std::endl;
    return 2;
  }

  Dune::GridPtr< GridType > gridPtr( argv[ 1 ] );

  GridView gridView = gridPtr->leafView();

  const unsigned int order   = atoi( argv[ 2 ] );
  const unsigned int level   = atoi( argv[ 3 ] );
  const unsigned int problem = atoi( argv[ 4 ] );

  Space space( gridView, order );
  const Space::DofMapper &dofMapper = space.dofMapper();
  std::vector< double > dofs( dofMapper.size() );

  Interpolation interpolation( space );
  interpolation( Function(problem), dofs );

  Dune::SubsamplingVTKWriter< GridView > vtkWriter( gridView, level );
  vtkWriter.addVertexData( new VTKFunction( space, dofs ) );
  vtkWriter.write( name );

  Dune::GrapeDataDisplay< GridType > grape( gridView );
  Dune::GrapeFunctionWrapper< Space > grapeFunction( space, dofs );
  grape.addData( grapeFunction.interface() );
  grape.display();
}
