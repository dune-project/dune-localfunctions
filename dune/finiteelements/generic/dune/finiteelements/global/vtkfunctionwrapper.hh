// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_VTKFUNCTIONWRAPPER_HH
#define DUNE_VTKFUNCTIONWRAPPER_HH

#include <dune/common/field.hh>

#include <dune/grid/io/file/vtk/vtkwriter.hh>

namespace Dune
{

  template< class DFS >
  class VTKFunctionWrapper
    : public VTKFunction< typename DFS::GridView::Grid >
  {
    typedef VTKFunctionWrapper< DFS > This;
    typedef VTKFunction< typename DFS::GridView::Grid > Base;

  public:
    typedef DFS DiscreteFunctionSpace;

    typedef typename DiscreteFunctionSpace::DomainField DomainField;
    typedef typename DiscreteFunctionSpace::DomainVector DomainVector;

    typedef typename DiscreteFunctionSpace::RangeField RangeField;
    typedef typename DiscreteFunctionSpace::RangeVector RangeVector;

    typedef typename DiscreteFunctionSpace::Basis Basis;
    typedef typename DiscreteFunctionSpace::GridView GridView;
    typedef typename GridView::template Codim< 0 >::Entity Entity;

    VTKFunctionWrapper ( const DiscreteFunctionSpace &dfSpace, const std::vector< RangeField > &dofs )
      : dfSpace_( dfSpace ),
        dofs_( dofs ),
        name_("VTKFunctionWrapper")
    {}
    VTKFunctionWrapper ( const std::string &name, const DiscreteFunctionSpace &dfSpace, const std::vector< RangeField > &dofs )
      : dfSpace_( dfSpace ),
        dofs_( dofs ),
        name_(name)
    {}

    virtual int ncomps () const
    {
      return DiscreteFunctionSpace::dimRange;
    }

    virtual double evaluate ( int comp, const Entity &entity, const DomainVector &x ) const
    {
      static std::vector< unsigned int > localIndices;
      dfSpace_.dofMapper().map( entity, localIndices );

      static std::vector< RangeVector > basisValues;
      const Basis &basis = dfSpace_.basis( entity );
      basisValues.resize( basis.size() );
      basis.evaluate( x, basisValues );

      RangeField y( 0 );
      for( unsigned int i = 0; i < basis.size(); ++i )
        y += dofs_[ localIndices[ i ] ] * basisValues[ i ][ comp ];
      return field_cast< double >( y );
    }

    virtual std::string name () const
    {
      return name_;
    }

  private:
    const DiscreteFunctionSpace &dfSpace_;
    const std::vector< RangeField > &dofs_;
    const std::string name_;
  };

}

#endif // #ifndef DUNE_VTKFUNCTIONWRAPPER_HH
