// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LAGRANGEPOINTS_HH
#define DUNE_LAGRANGEPOINTS_HH

#include <vector>

#include <dune/common/fvector.hh>

#include <dune/grid/genericgeometry/topologytypes.hh>
#include <dune/grid/genericgeometry/subtopologies.hh>

#include <dune/finiteelements/common/localcoefficients.hh>

namespace Dune
{

  // Internal Forward Declarations
  // -----------------------------

  template< class Topology, class F >
  class LagrangePoints;



  // LagrangePoint
  // -------------

  template< class F, unsigned int dim >
  class LagrangePoint
  {
    typedef LagrangePoint< F, dim > This;

    template< class, class >
    friend class LagrangePointsImpl;

  public:
    static const int dimension = dim;

    typedef F Field;

    typedef FieldVector< Field, dimension > Vector;

  private:
    Vector point_;
    LocalKey localKey_;

  public:
    const Vector &point () const
    {
      return point_;
    }

    const LocalKey &localKey () const
    {
      return localKey_;
    }
  };



  // LagrangePointsImpl
  // -------------....-

  template< class Topology, class F >
  class LagrangePointsImpl;

  template< class F >
  class LagrangePointsImpl< GenericGeometry::Point, F >
  {
    typedef LagrangePointsImpl< GenericGeometry::Point, F > This;

    typedef GenericGeometry::Point Topology;

  public:
    typedef F Field;

    static const unsigned int dimension = Topology::dimension;

  private:
    friend class LagrangePoints< Topology, Field >;
    friend class LagrangePointsImpl< GenericGeometry::Prism< Topology >, Field >;
    friend class LagrangePointsImpl< GenericGeometry::Pyramid< Topology >, Field >;

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

    static unsigned int size ( const unsigned int order )
    {
      return 1;
    }
  };

  template< class BaseTopology, class F >
  class LagrangePointsImpl< GenericGeometry::Prism< BaseTopology >, F >
  {
    typedef LagrangePointsImpl< GenericGeometry::Prism< BaseTopology >, F > This;

    typedef GenericGeometry::Prism< BaseTopology > Topology;

  public:
    typedef F Field;

    static const unsigned int dimension = Topology::dimension;

  private:
    friend class LagrangePoints< Topology, Field >;
    friend class LagrangePointsImpl< GenericGeometry::Prism< Topology >, Field >;
    friend class LagrangePointsImpl< GenericGeometry::Pyramid< Topology >, Field >;

    typedef LagrangePointsImpl< BaseTopology, Field > BaseImpl;

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
            key = LocalKey( key.subentity(), codim, key.index() );
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
          key = LocalKey( key.subentity() + numBaseN, codim, key.index() );
          points[ j + n ].point_ = points[ j ].point_;
          points[ j + n ].point_[ dimension-1 ] = Field( 1 );
          points[ j + n ].localKey_ = LocalKey( key.subentity() + numBaseM, codim, key.index() );
          ++count[ key.subentity() + numBaseM ];
        }
        size += 2*n;
      }
      return size;
    }

    static unsigned int size ( const unsigned int order )
    {
      return BaseImpl::size( order ) * (order+1);
    }
  };

  template< class BaseTopology, class F >
  class LagrangePointsImpl< GenericGeometry::Pyramid< BaseTopology >, F >
  {
    typedef LagrangePointsImpl< GenericGeometry::Pyramid< BaseTopology >, F > This;

    typedef GenericGeometry::Pyramid< BaseTopology > Topology;

  public:
    typedef F Field;

    static const unsigned int dimension = Topology::dimension;

  private:
    friend class LagrangePoints< Topology, Field >;
    friend class LagrangePointsImpl< GenericGeometry::Prism< Topology >, Field >;
    friend class LagrangePointsImpl< GenericGeometry::Pyramid< Topology >, Field >;

    typedef LagrangePointsImpl< BaseTopology, Field > BaseImpl;

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
          key = LocalKey( key.subentity(), codim, key.index() );
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
            key = LocalKey( key.subentity()+numBaseM, codim, key.index() );
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

    static unsigned int size ( const unsigned int order )
    {
      unsigned int size = BaseImpl::size( order );
      for( unsigned int i = 1; i <= order; ++i )
        size += BaseImpl::size( order - i );
      return size;
    }
  };



  // LagnragePoints
  // --------------
  template< class Topology, class F >
  class LagrangePoints
  {
    typedef LagrangePoints< Topology, F > This;
    typedef LagrangePointsImpl< Topology, F > Impl;

  public:
    typedef F Field;

    static const unsigned int dimension = Topology::dimension;

    typedef Dune::LagrangePoint< Field, dimension > LagrangePoint;

    typedef typename std::vector< LagrangePoint >::const_iterator iterator;

  private:
    template< int codim >
    struct Init;

    unsigned int order_;
    std::vector< LagrangePoint > points_;

  public:
    LagrangePoints ( const unsigned int order )
      : order_( order ),
        points_( Impl::size( order ) )
    {
      LagrangePoint *p = &(points_[ 0 ]);
      GenericGeometry::ForLoop< Init, 0, dimension >::apply( order, p );
    }

    iterator begin () const
    {
      return points_.begin();
    }

    iterator end () const
    {
      return points_.end();
    }

    unsigned int order () const
    {
      return order_;
    }

    size_t size () const
    {
      return points_.size();
    }
  };



  // LagrangePoints::Init

  template< class Topology, class F >
  template< int dim >
  struct LagrangePoints< Topology, F >::Init
  {
    static void apply ( const unsigned int order, LagrangePoint *&p )
    {
      const unsigned int size = GenericGeometry::Size< Topology, dimension-dim >::value;
      unsigned int count[ size ];
      for( unsigned int i = 0; i < size; ++i )
        count[ i ] = 0;
      p += Impl::template setup< dimension-dim, dimension >( order, count, p );
    }
  };

}

#endif
