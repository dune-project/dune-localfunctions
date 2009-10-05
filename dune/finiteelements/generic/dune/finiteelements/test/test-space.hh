// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <dune/grid/genericgeometry/topologytypes.hh>
#if HAVE_DUNE_PSG
#include <dune/grid/io/file/dgfparser/dgfpsggridtype.hh>
#endif
#include <dune/grid/io/file/dgfparser/dgfgridtype.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/grid/io/visual/grapedatadisplay.hh>

#include <dune/finiteelements/global/vtkfunctionwrapper.hh>
#include <dune/finiteelements/global/grapefunctionwrapper.hh>

const unsigned int dimension = GridType::dimension;

typedef GridType::LeafGridView GridView;

typedef GridView::Codim< 0 >::Entity Entity;
typedef GridView::Codim< 0 >::Iterator Iterator;

struct Function
{
  typedef Dune::FieldVector< double, dimension > DomainVector;
  typedef Dune::FieldVector< double, 1 > RangeVector;

  Function(unsigned int p) : problem_(p) {}

  RangeVector operator() ( const DomainVector &x ) const
  {
    if (problem_ == 1)
      if (x[0]<0.5)
        return exp( -x.two_norm() );
      else
        return (1.+sin( x.two_norm()*2.*M_PI ))*0.2;
    else if (problem_ == 2)
      if (x[0]<x[1])
        return exp( -x.two_norm() );
      else
        return (1.+sin( x.two_norm()*2.*M_PI ))*0.2;
    else if (problem_ == 3)
      if (x.two_norm()<0.6)
        return exp( -x.two_norm() );
      else
        return (1.+sin( x.two_norm()*2.*M_PI ))*0.2;
    else
      return exp( -x.two_norm2() );
  }
  unsigned int problem_;
};
