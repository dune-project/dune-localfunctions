// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRAPEFUNCTIONWRAPPER_HH
#define DUNE_GRAPEFUNCTIONWRAPPER_HH

#include <dune/finiteelements/common/field.hh>

#if HAVE_GRAPE
#include <dune/grid/io/visual/grapedatadisplay.hh>
#endif // #if HAVE_GRAPE

namespace Dune
{

#if not HAVE_GRAPE
  template< class GV, int dimR, int polOrd >
  struct GrapeFunction
  {};
#endif



  // GrapeFunctionWrapper
  // -------------------.

  template< class DFS >
  class GrapeFunctionWrapper
    : public GrapeFunction< typename DFS::GridView, 1, 4 >
  {
    typedef GrapeFunctionWrapper< DFS > This;
    typedef GrapeFunction< typename DFS::GridView, 1, 4 > Base;

  public:
    typedef DFS DiscreteFunctionSpace;

    typedef typename DiscreteFunctionSpace::DomainField DomainField;
    typedef typename DiscreteFunctionSpace::DomainVector DomainVector;

    typedef typename DiscreteFunctionSpace::RangeField RangeField;
    typedef typename DiscreteFunctionSpace::RangeVector RangeVector;

    typedef typename DiscreteFunctionSpace::Basis Basis;
    typedef typename DiscreteFunctionSpace::GridView GridView;
    typedef typename GridView::template Codim< 0 >::Entity Entity;

    GrapeFunctionWrapper ( const DiscreteFunctionSpace &dfSpace, const std::vector< double > &dofs )
      : dfSpace_( dfSpace ),
        dofs_( dofs )
    {}

    virtual void evaluate ( const Entity &entity, const DomainVector &x, RangeVector &y ) const
    {
      static std::vector< unsigned int > localIndices;
      dfSpace_.dofMapper().map( entity, localIndices );

      static std::vector< RangeVector > basisValues;
      const Basis &basis = dfSpace_.basis( entity );
      basisValues.resize( basis.size() );
      basis.evaluate( x, basisValues );

      y = RangeVector( Zero< RangeField >() );
      for( unsigned int i = 0; i < basis.size(); ++i )
        y += dofs_[ localIndices[ i ] ] * basisValues[ i ];
    }

    virtual const GridView &gridView () const
    {
      return dfSpace_.gridView();
    }

    virtual std::string name () const
    {
      return "GrapeFunctionWrapper";
    }

    const Base &interface () const
    {
      return *this;
    }

  private:
    const DiscreteFunctionSpace &dfSpace_;
    const std::vector< RangeField > &dofs_;
  };

}

#endif // #ifndef DUNE_GRAPEFUNCTIONWRAPPER_HH
