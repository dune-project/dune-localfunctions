// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LAGRANGEPOINTS_HH
#define DUNE_LAGRANGEPOINTS_HH

#include <vector>

#include <dune/common/fvector.hh>

#include <dune/grid/genericgeometry/topologytypes.hh>

namespace Dune
{

  // Internal Forward Declarations
  // -----------------------------

  template< class Topology, class F >
  class LagrangePoints;



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

    template< unsigned int dim >
    static FieldVector< Field, dim > *
    setup ( const unsigned int order, FieldVector< Field, dim > *points )
    {
      *(points++) = Field( 0 );
      return points;
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

    template< unsigned int dim >
    static FieldVector< Field, dim > *
    setup ( const unsigned int order, FieldVector< Field, dim > *points )
    {
      for( unsigned int i = 0; i <= order; ++i )
      {
        FieldVector< Field, dim > *const end
          = BaseImpl::template setup< dim >( order, points );
        for( ; points != end; ++points )
          (*points)[ dimension-1 ] = Field( i ) / Field( order );
      }
      return points;
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

    template< unsigned int dim >
    static FieldVector< Field, dim > *
    setup ( const unsigned int order, FieldVector< Field, dim > *points )
    {
      if( order > 0 )
      {
        for( unsigned int i = 0; i <= order; ++i )
        {
          FieldVector< Field, dim > *const end
            = BaseImpl::template setup< dim >( order - i, points );
          for( ; points != end; ++points )
          {
            for( unsigned int j = 0; j < dimension-1; ++j )
              (*points)[ j ] *= Field( order - i ) / Field( order );
            (*points)[ dimension-1 ] = Field( i ) / Field( order );
          }
        }
      }
      else
        *(points++) = Field( 0 );
      return points;
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

    typedef FieldVector< Field, dimension > Vector;

    typedef typename std::vector< Vector >::const_iterator iterator;

  private:
    unsigned int order_;
    std::vector< Vector > points_;

  public:
    LagrangePoints ( unsigned int order )
      : order_( order ),
        points_( Impl::size( order ) )
    {
      Impl::template setup< dimension >( order, &(points_[ 0 ]) );
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

}

#endif
