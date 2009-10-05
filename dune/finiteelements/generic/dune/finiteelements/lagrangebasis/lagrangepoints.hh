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

  template< class F, unsigned int dim >
  class LagrangePoint;

  template< class F, unsigned int dim >
  class LagrangePoints;

  template< class F, unsigned int dim >
  class LagrangePointsCreator;



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

    const Vector &point () const
    {
      return point_;
    }

    const LocalKey &localKey () const
    {
      return localKey_;
    }

  private:
    Vector point_;
    LocalKey localKey_;
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

    friend class LagrangePointsCreator< F, Topology::dimension >;
    friend class LagrangePointsImpl< GenericGeometry::Prism< Topology >, F >;
    friend class LagrangePointsImpl< GenericGeometry::Pyramid< Topology >, F >;

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
  class LagrangePointsImpl< GenericGeometry::Prism< BaseTopology >, F >
  {
    typedef LagrangePointsImpl< GenericGeometry::Prism< BaseTopology >, F > This;

    typedef GenericGeometry::Prism< BaseTopology > Topology;

    friend class LagrangePointsCreator< F, Topology::dimension >;
    friend class LagrangePointsImpl< GenericGeometry::Prism< Topology >, F >;
    friend class LagrangePointsImpl< GenericGeometry::Pyramid< Topology >, F >;

    typedef LagrangePointsImpl< BaseTopology, F > BaseImpl;

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
  class LagrangePointsImpl< GenericGeometry::Pyramid< BaseTopology >, F >
  {
    typedef LagrangePointsImpl< GenericGeometry::Pyramid< BaseTopology >, F > This;

    typedef GenericGeometry::Pyramid< BaseTopology > Topology;

    friend class LagrangePointsCreator< F, Topology::dimension >;
    friend class LagrangePointsImpl< GenericGeometry::Prism< Topology >, F >;
    friend class LagrangePointsImpl< GenericGeometry::Pyramid< Topology >, F >;

    typedef LagrangePointsImpl< BaseTopology, F > BaseImpl;

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
  class LagrangePoints
    : public LocalCoefficientsInterface< LagrangePoints< F, dim > >
  {
    typedef LagrangePoints< F, dim > This;

    friend class LagrangePointsCreator< F, dim >;

  public:
    typedef F Field;

    static const unsigned int dimension = dim;

    typedef Dune::LagrangePoint< Field, dimension > LagrangePoint;

    typedef typename std::vector< LagrangePoint >::const_iterator iterator;

  private:
    LagrangePoints ( const unsigned int order, const unsigned int size )
      : order_( order ),
        points_( size )
    {}

  public:
    const LagrangePoint &operator[] ( const unsigned int i ) const
    {
      assert( i < size() );
      return points_[ i ];
    }

    iterator begin () const
    {
      return points_.begin();
    }

    iterator end () const
    {
      return points_.end();
    }

    const LocalKey &localKey ( const unsigned int i ) const
    {
      return (*this)[ i ].localKey();
    }

    unsigned int order () const
    {
      return order_;
    }

    unsigned int size () const
    {
      return points_.size();
    }

  private:
    unsigned int order_;
    std::vector< LagrangePoint > points_;
  };



  // LagrangePointsCreator
  // ---------------------

  template< class F, unsigned int dim >
  class LagrangePointsCreator
  {
    template< class T >
    struct Topology;

  public:
    static const unsigned int dimension = dim;

    typedef Dune::LagrangePoints< F, dimension > LagrangePoints;
    typedef LagrangePoints LocalCoefficients;

    typedef unsigned int Key;

    template< class T >
    static const LagrangePoints &lagrangePoints ( const Key &order );

    template< class T >
    static const LocalCoefficients &localCoefficients ( const Key &order )
    {
      return lagrangePoints< T >( order );
    }

    static void release ( const LagrangePoints &lagrangePoints )
    {
      delete &lagrangePoints;
    }
  };


  template< class F, unsigned int dim >
  template< class T >
  struct LagrangePointsCreator< F, dim >::Topology
  {
    typedef LagrangePointsImpl< T, F > Impl;
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
  inline const typename LagrangePointsCreator< F, dim >::LagrangePoints &
  LagrangePointsCreator< F, dim >::lagrangePoints ( const Key &order )
  {
    typedef Dune::LagrangePoint< F, dimension > LagrangePoint;
    typedef typename Topology< T >::Impl Impl;

    LagrangePoints *lagrangePoints = new LagrangePoints( order, Impl::size( order ) );
    LagrangePoint *p = &(lagrangePoints->points_[ 0 ]);
    GenericGeometry::ForLoop< Topology< T >::template Init, 0, dimension >::apply( order, p );
    return *lagrangePoints;
  }

}

#endif
