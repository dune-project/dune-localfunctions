// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_LAGRANGE_EQUIDISTANTPOINTS_HH
#define DUNE_LOCALFUNCTIONS_LAGRANGE_EQUIDISTANTPOINTS_HH

#include <cstddef>

#include <algorithm>
#include <vector>

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>

#include <dune/localfunctions/lagrange/emptypoints.hh>
#include <dune/localfunctions/utility/field.hh>

namespace Dune
{

  // numLagrangePoints
  // -----------------

  inline std::size_t numLagrangePoints ( const GeometryType& gt, std::size_t order )
  {
    const int dim = gt.dim();
    if( dim > 0 )
    {
      const GeometryType baseGeometryType = Impl::getBase( gt );
      if( gt.isConical() )
      {
        std::size_t size = 0;
        for( unsigned int o = 0; o <= order; ++o )
          size += numLagrangePoints( baseGeometryType, o );
        return size;
      }
      else
        return numLagrangePoints( baseGeometryType, order ) * (order+1);
    }
    else
      return 1;
  }

  [[deprecated("Use numLagrangePoints(const GeometryType& gt, std::size_t order ) instead.")]]
  inline std::size_t numLagrangePoints (  unsigned int topologyId, unsigned int dim, std::size_t order )
  {
    return numLagrangePoints ( GeometryType(topologyId, dim), order);
  }



  // equidistantLagrangePoints
  // -------------------------

  template< class ct, unsigned int cdim >
  inline static unsigned int equidistantLagrangePoints ( const GeometryType& gt, unsigned int codim, std::size_t order, unsigned int *count, LagrangePoint< ct, cdim > *points )
  {
    const unsigned int dim = gt.dim();
    assert( (0 <= codim) && (codim <= dim) && (dim <= cdim) );

    if( dim > 0 )
    {
      const GeometryType baseGeometryType = Impl::getBase( gt );
      const unsigned int numBaseN = (codim < dim ? Geo::Impl::size( baseGeometryType.id(), baseGeometryType.dim(), codim ) : 0);
      const unsigned int numBaseM = (codim > 0 ? Geo::Impl::size( baseGeometryType.id(), baseGeometryType.dim(), codim-1 ) : 0);

      if( gt.isPrismatic() )
      {
        unsigned int size = 0;
        if( codim < dim )
        {
          for( unsigned int i = 1; i < order; ++i )
          {
            const unsigned int n = equidistantLagrangePoints( baseGeometryType, codim, order, count, points );
            for( unsigned int j = 0; j < n; ++j )
            {
              LocalKey &key = points->localKey_;
              key = LocalKey( key.subEntity(), codim, key.index() );
              points->point_[ dim-1 ] = ct( i ) / ct( order );
              ++points;
            }
            size += n;
          }
        }

        if( codim > 0 )
        {
          const unsigned int n = equidistantLagrangePoints( baseGeometryType, codim-1, order, count+numBaseN, points );
          for( unsigned int j = 0; j < n; ++j )
          {
            LocalKey &key = points[ j ].localKey_;
            key = LocalKey( key.subEntity() + numBaseN, codim, key.index() );

            points[ j + n ].point_ = points[ j ].point_;
            points[ j + n ].point_[ dim-1 ] = ct( 1 );
            points[ j + n ].localKey_ = LocalKey( key.subEntity() + numBaseM, codim, key.index() );
            ++count[ key.subEntity() + numBaseM ];
          }
          size += 2*n;
        }

        return size;
      }
      else
      {
        unsigned int size = (codim > 0 ? equidistantLagrangePoints( baseGeometryType, codim-1, order, count, points ) : 0);
        LagrangePoint< ct, cdim > *const end = points + size;
        for( ; points != end; ++points )
          points->localKey_ = LocalKey( points->localKey_.subEntity(), codim, points->localKey_.index() );

        if( codim < dim )
        {
          for( unsigned int i = order-1; i > 0; --i )
          {
            const unsigned int n = equidistantLagrangePoints( baseGeometryType, codim, i, count+numBaseM, points );
            LagrangePoint< ct, cdim > *const end = points + n;
            for( ; points != end; ++points )
            {
              points->localKey_ = LocalKey( points->localKey_.subEntity()+numBaseM, codim, points->localKey_.index() );
              for( unsigned int j = 0; j < dim-1; ++j )
                points->point_[ j ] *= ct( i ) / ct( order );
              points->point_[ dim-1 ] = ct( order - i ) / ct( order );
            }
            size += n;
          }
        }
        else
        {
          points->localKey_ = LocalKey( numBaseM, dim, count[ numBaseM ]++ );
          points->point_ = 0;
          points->point_[ dim-1 ] = ct( 1 );
          ++size;
        }

        return size;
      }
    }
    else
    {
      points->localKey_ = LocalKey( 0, 0, count[ 0 ]++ );
      points->point_ = 0;
      return 1;
    }
  }

  template< class ct, unsigned int cdim >
  [[deprecated("Use equidistantLagrangePoints ( GeometryType gt,  ... ) instead.")]]
  inline static unsigned int equidistantLagrangePoints ( unsigned int topologyId, unsigned int dim, unsigned int codim, std::size_t order, unsigned int *count, LagrangePoint< ct, cdim > *points )
  {
    return equidistantLagrangePoints ( GeometryType(topologyId, dim), codim, order, *count, *points );
  }



  // EquidistantPointSet
  // -------------------

  template< class F, unsigned int dim >
  class EquidistantPointSet
    : public EmptyPointSet< F, dim >
  {
    typedef EmptyPointSet< F, dim > Base;

  public:
    static const unsigned int dimension = dim;

    using Base::order;

    EquidistantPointSet ( std::size_t order ) : Base( order ) {}

    void build ( GeometryType gt )
    {
      assert( gt.dim() == dimension );
      points_.resize( numLagrangePoints( gt, order() ) );

      typename Base::LagrangePoint *p = points_.data();
      std::vector< unsigned int > count;
      for( unsigned int mydim = 0; mydim <= dimension; ++mydim )
      {
        count.resize( Geo::Impl::size( gt.id(), dimension, dimension-mydim ) );
        std::fill( count.begin(), count.end(), 0u );
        p += equidistantLagrangePoints( gt, dimension-mydim, order(), count.data(), p );
      }
      const auto &refElement = referenceElement<F,dimension>(gt);
      F weight = refElement.volume()/F(double(points_.size()));
      for (auto &p : points_)
        p.weight_ = weight;
    }

    template< GeometryType::Id geometryId >
    bool build ()
    {
      build( GeometryType( geometryId ) );
      return true;
    }

    bool buildCube ()
    {
      return build< GeometryTypes::cube(dim) > ();
    }

    static bool supports ( GeometryType, std::size_t /*order*/ ) { return true; }
    template< GeometryType::Id geometryId>
    static bool supports ( std::size_t order ) {
      return supports( GeometryType( geometryId ), order );
    }

  private:
    using Base::points_;
  };

} // namespace Dune

#endif // #ifndef DUNE_LOCALFUNCTIONS_LAGRANGE_EQUIDISTANTPOINTS_HH
