// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <dune/grid/genericgeometry/topologytypes.hh>
#include <dune/grid/io/file/dgfparser/dgfgridtype.hh>

#include <dune/finiteelements/orthonormalbasis/dgorthonormalbasis.hh>
#include <dune/finiteelements/global/dofmapper.hh>
#include <dune/finiteelements/global/interpolation.hh>

const unsigned int dimension = GridType::dimension;

typedef double StorageField;
typedef Dune::AlgLib::MultiPrecision< 512 > ComputeField;
typedef Dune::DGOrthonormalBasisProvider< dimension, StorageField, ComputeField > BasisProvider;

typedef GridType::LeafGridView GridView;
typedef Dune::DofMapper< GridView::IndexSet, BasisProvider > DofMapper;

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

  DofMapper dofMapper( gridView.indexSet(), order );

  const Iterator end = gridView.end< 0 >();
  for( Iterator it = gridView.begin< 0 >(); it != end; ++it )
  {
    const Entity &entity = *it;
    const unsigned int topologyId = Dune::GenericGeometry::topologyId( entity.type() );

    const BasisProvider::Basis &basis = BasisProvider::basis( topologyId, order );
  }
}
