// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_EQUIDISTANTPOINTS_HH
#define DUNE_EQUIDISTANTPOINTS_HH

#include <vector>

#include <dune/common/fvector.hh>
#include <dune/finiteelements/common/field.hh>
#include <dune/common/forloop.hh>
#include <dune/finiteelements/generic/topologyfactory.hh>
#include <dune/grid/genericgeometry/topologytypes.hh>
#include <dune/grid/genericgeometry/subtopologies.hh>

#include <dune/finiteelements/lagrangebasis/emptypoints.hh>

namespace Dune
{

  // Internal Forward Declarations
  // -----------------------------

  template< class F, unsigned int dim >
  class EquidistantPointSet;

  // EquidistantPointSetImpl
  // ----------------------------

  template< class Topology, class F >
  class EquidistantPointSetImpl;

  template< class F >
  class EquidistantPointSetImpl< GenericGeometry::Point, F >
  {
    typedef EquidistantPointSetImpl< GenericGeometry::Point, F > This;

    typedef GenericGeometry::Point Topology;

    friend class EquidistantPointSet< F, Topology::dimension >;
    friend class EquidistantPointSetImpl< GenericGeometry::Prism< Topology >, F >;
    friend class EquidistantPointSetImpl< GenericGeometry::Pyramid< Topology >, F >;

  public:
    typedef F Field;

    static const unsigned int dimension = Topology::dimension;

    static unsigned int size ( const unsigned int order )
    {
      return 1;
    }

  private:
    template< unsigned int codim, unsigned int dim >
    static unsigned int setup ( const unsigned int order,
                                unsigned int *count,
                                LagrangePoint< Field, dim > *points )
    {
      assert( codim == 0 );
      points->localKey_ = LocalKey( 0, 0, count[ 0 ]++ );
      points->point_ = Field( 0 );
      return 1;
    }
  };

  template< class BaseTopology, class F >
  class EquidistantPointSetImpl< GenericGeometry::Prism< BaseTopology >, F >
  {
    typedef EquidistantPointSetImpl< GenericGeometry::Prism< BaseTopology >, F > This;

    typedef GenericGeometry::Prism< BaseTopology > Topology;

    friend class EquidistantPointSet< F, Topology::dimension >;
    friend class EquidistantPointSetImpl< GenericGeometry::Prism< Topology >, F >;
    friend class EquidistantPointSetImpl< GenericGeometry::Pyramid< Topology >, F >;

    typedef EquidistantPointSetImpl< BaseTopology, F > BaseImpl;

  public:
    typedef F Field;

    static const unsigned int dimension = Topology::dimension;

    static unsigned int size ( const unsigned int order )
    {
      return BaseImpl::size( order ) * (order+1);
    }

    // private:
    template< unsigned int codim, unsigned int dim >
    static unsigned int setup ( const unsigned int order,
                                unsigned int *count,
                                LagrangePoint< Field, dim > *points )
    {
      unsigned int size = 0;
      unsigned int numBaseN = 0;

      if( codim < dimension )
      {
        const unsigned int vcodim = (codim < dimension ? codim : dimension-1);
        numBaseN = GenericGeometry::Size< BaseTopology, vcodim >::value;
        for( unsigned int i = 1; i < order; ++i )
        {
          const unsigned int n = BaseImpl::template setup< vcodim, dim >( order, count, points );
          for( unsigned int j = 0; j < n; ++j )
          {
            LocalKey &key = points->localKey_;
            key = LocalKey( key.subEntity(), codim, key.index() );
            points->point_[ dimension-1 ] = Field( i ) / Field( order );
            ++points;
          }
          size += n;
        }
      }

      if( codim > 0 )
      {
        const unsigned int vcodim = (codim > 0 ? codim : 1);
        const unsigned int numBaseM = GenericGeometry::Size< BaseTopology, vcodim-1 >::value;
        const unsigned int n = BaseImpl::template setup< vcodim-1, dim >( order, count+numBaseN, points );
        for( unsigned int j = 0; j < n; ++j )
        {
          LocalKey &key = points[ j ].localKey_;
          key = LocalKey( key.subEntity() + numBaseN, codim, key.index() );
          points[ j + n ].point_ = points[ j ].point_;
          points[ j + n ].point_[ dimension-1 ] = Field( 1 );
          points[ j + n ].localKey_ = LocalKey( key.subEntity() + numBaseM, codim, key.index() );
          ++count[ key.subEntity() + numBaseM ];
        }
        size += 2*n;
      }
      return size;
    }
  };

  template< class BaseTopology, class F >
  class EquidistantPointSetImpl< GenericGeometry::Pyramid< BaseTopology >, F >
  {
    typedef EquidistantPointSetImpl< GenericGeometry::Pyramid< BaseTopology >, F > This;

    typedef GenericGeometry::Pyramid< BaseTopology > Topology;

    friend class EquidistantPointSet< F, Topology::dimension >;
    friend class EquidistantPointSetImpl< GenericGeometry::Prism< Topology >, F >;
    friend class EquidistantPointSetImpl< GenericGeometry::Pyramid< Topology >, F >;

    typedef EquidistantPointSetImpl< BaseTopology, F > BaseImpl;

  public:
    typedef F Field;

    static const unsigned int dimension = Topology::dimension;

    static unsigned int size ( const unsigned int order )
    {
      unsigned int size = BaseImpl::size( order );
      for( unsigned int i = 1; i <= order; ++i )
        size += BaseImpl::size( order - i );
      return size;
    }

    // private:
    template< unsigned int codim, unsigned int dim >
    static unsigned int setup ( const unsigned int order,
                                unsigned int *count,
                                LagrangePoint< Field, dim > *points )
    {
      unsigned int size = 0;
      unsigned int numBaseM = 0;

      if( codim > 0 )
      {
        const unsigned int vcodim = (codim > 0 ? codim : 1);
        numBaseM = GenericGeometry::Size< BaseTopology, vcodim-1 >::value;
        size = BaseImpl::template setup< vcodim-1, dim >( order, count, points );
        LagrangePoint< Field, dim > *const end = points + size;
        for( ; points != end; ++points )
        {
          LocalKey &key = points->localKey_;
          key = LocalKey( key.subEntity(), codim, key.index() );
        }
      }

      if( codim < dimension )
      {
        const unsigned int vcodim = (codim < dimension ? codim : dimension-1);
        for( unsigned int i = order-1; i > 0; --i )
        {
          const unsigned int n = BaseImpl::template setup< vcodim, dim >( i, count+numBaseM, points );
          LagrangePoint< Field, dim > *const end = points + n;
          for( ; points != end; ++points )
          {
            LocalKey &key = points->localKey_;
            key = LocalKey( key.subEntity()+numBaseM, codim, key.index() );
            for( unsigned int j = 0; j < dimension-1; ++j )
              points->point_[ j ] *= Field( i ) / Field( order );
            points->point_[ dimension-1 ] = Field( order - i ) / Field( order );
          }
          size += n;
        }
      }
      else
      {
        points->localKey_ = LocalKey( numBaseM, dimension, count[ numBaseM ]++ );
        points->point_ = Field( 0 );
        points->point_[ dimension-1 ] = 1;
        ++size;
      }

      return size;
    }
  };



  // LagnragePoints
  // --------------

  template< class F, unsigned int dim >
  class EquidistantPointSet : public EmptyPointSet<F,dim>
  {
    template< class T >
    struct Topology;
    typedef EmptyPointSet<F,dim> Base;
  public:
    static const unsigned int dimension = dim;

    EquidistantPointSet( unsigned int order )
      : Base(order)
    {}

    template< class T >
    bool build ( );

    template< class T >
    static bool supports ( unsigned int order )
    {
      return true;
    }
  private:
    using Base::points_;
  };

  template< class F, unsigned int dim >
  template< class T >
  struct EquidistantPointSet< F, dim >::Topology
  {
    typedef EquidistantPointSetImpl< T, F > Impl;
    typedef Dune::LagrangePoint< F, dim > LagrangePoint;

    template< int pdim >
    struct Init
    {
      static void apply ( const unsigned int order, LagrangePoint *&p )
      {
        const unsigned int size = GenericGeometry::Size< T, dimension-pdim >::value;
        unsigned int count[ size ];
        for( unsigned int i = 0; i < size; ++i )
          count[ i ] = 0;
        p += Impl::template setup< dimension-pdim, dimension >( order, count, p );
      }
    };
  };

  template< class F, unsigned int dim >
  template< class T >
  inline bool EquidistantPointSet< F, dim >::build ( )
  {
    unsigned int order = Base::order();
    typedef Dune::LagrangePoint< F, dimension > LagrangePoint;
    typedef typename Topology< T >::Impl Impl;
    points_.resize( Impl::size( order ) );
    LagrangePoint *p = &(points_[ 0 ]);
    ForLoop< Topology< T >::template Init, 0, dimension >::apply( order, p );
    return true;
  }
}

#endif
