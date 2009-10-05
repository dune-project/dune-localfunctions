// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <dune/grid/genericgeometry/topologytypes.hh>
#include <dune/grid/io/file/dgfparser/dgfgridtype.hh>

#include <dune/finiteelements/orthonormalbasis/dgorthonormalbasis.hh>
#include <dune/finiteelements/dofmapper.hh>

const unsigned int dimension = GridType::dimension;

typedef double StorageField;
typedef Dune::AlgLib::MultiPrecision< 512 > ComputeField;
typedef Dune::DGOrthonormalBasisProvider< dimension, StorageField, ComputeField > BasisProvider;

typedef Dune::DofMapper< GridType::LeafIndexSet, BasisProvider > DofMapper;

typedef GridType::Codim< 0 >::Entity Entity;
typedef GridType::Codim< 0 >::LeafIterator LeafIterator;

int main ( int argc, char **argv )
{
  if( argc < 3 )
  {
    std::cerr << "Usage: " << argv[ 0 ] << " <dgf-file> <order>" << std::endl;
    return 2;
  }

  Dune::GridPtr< GridType > gridPtr( argv[ 1 ] );
  GridType &grid = *gridPtr;

  const unsigned int order = atoi( argv[ 2 ] );

  DofMapper dofMapper( grid.leafIndexSet(), order );

  const LeafIterator end = grid.leafend< 0 >();
  for( LeafIterator it = grid.leafbegin< 0 >(); it != end; ++it )
  {
    const Entity &entity = *it;
    const unsigned int topologyId = Dune::GenericGeometry::topologyId( entity.type() );

    const BasisProvider::Basis &basis = BasisProvider::basis( topologyId, order );
  }
}
