// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <dune/grid/genericgeometry/topologytypes.hh>
#include <dune/grid/io/file/dgfparser/dgfgridtype.hh>

#include <dune/finiteelements/lagrangebasis/space.hh>
#include <dune/finiteelements/global/vtkfunctionwrapper.hh>

const unsigned int dimension = GridType::dimension;

typedef GridType::LeafGridView GridView;
typedef Dune::LagrangeSpace< GridView, double > Space;
typedef Dune::VTKFunctionWrapper< Space > VTKFunction;

typedef GridView::Codim< 0 >::Entity Entity;
typedef GridView::Codim< 0 >::Iterator Iterator;

int main ( int argc, char **argv )
{
  if( argc < 3 )
  {
    std::cerr << "Usage: " << argv[ 0 ] << " <dgf-file> <order>" << std::endl;
    return 2;
  }

  Dune::GridPtr< GridType > gridPtr( argv[ 1 ] );

  GridView gridView = gridPtr->leafView();

  const unsigned int order = atoi( argv[ 2 ] );

  Space space( gridView, order );
  const Space::DofMapper &dofMapper = space.dofMapper();
  std::vector< double > dofs( dofMapper.size() );
  VTKFunction vtkFunction( space, dofs );

  const Iterator end = gridView.end< 0 >();
  for( Iterator it = gridView.begin< 0 >(); it != end; ++it )
  {
    const Entity &entity = *it;

    const Space::Basis &basis = space.basis( entity );

  }
}
