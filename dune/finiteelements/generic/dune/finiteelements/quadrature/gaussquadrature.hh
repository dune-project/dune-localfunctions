// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GAUSSQUADRATURE_HH
#define DUNE_GAUSSQUADRATURE_HH

#include <alglib/gqgengauss.h>

#include <dune/common/field.hh>

#include <dune/alglib/multiprecision.hh>
#include <dune/alglib/vector.hh>

#include <dune/finiteelements/quadrature/quadrature.hh>

namespace Dune
{

  namespace GenericGeometry
  {

    // GaussPoints
    // -----------

    template< class F >
    class GaussPoints;

    template< unsigned int precision >
    class GaussPoints< AlgLib::MultiPrecision< precision > >
    {
      typedef GaussPoints< AlgLib::MultiPrecision< precision > > This;

    public:
      typedef AlgLib::MultiPrecision< precision > Field;
      typedef AlgLib::Vector< Field > Vector;

      struct Iterator;

      explicit GaussPoints ( unsigned int n )
        : points_( n ),
          weights_( n )
      {
        Vector alpha( n );
        Vector beta( n );
        for( unsigned int i = 0; i < n; ++i )
        {
          alpha[ i ] = 0;
          beta[ i ] = Field( i*i ) / Field( 4*(i*i)-1 );
        }

        gqgengauss::generategaussquadrature< precision >( alpha, beta, Field( 2 ), n, points_, weights_ );

        const Field half = Field( 1 ) / Field( 2 );
        for( unsigned int i = 0; i < n; ++i )
        {
          points_[ i ] = (points_[ i ] + Field( 1 )) * half;
          weights_[ i ] *= half;
        }
      }

      Iterator begin () const
      {
        return Iterator( *this, 0 );
      }

      Iterator end () const
      {
        return Iterator( *this, size() );
      }

      Field point ( const unsigned int i ) const
      {
        assert( i < points_.size() );
        return points_[ i ];
      }

      Field weight ( const unsigned int i ) const
      {
        assert( i < weights_.size() );
        return weights_[ i ];
      }

      unsigned int size () const
      {
        assert( points_.size() == weights_.size() );
        return points_.size();
      }

    private:
      Vector points_;
      Vector weights_;
    };



    // GaussPoints::Iterator
    // ---------------------

    template< unsigned int precision >
    struct GaussPoints< AlgLib::MultiPrecision< precision > >::Iterator
    {
      typedef AlgLib::MultiPrecision< precision > Field;
      typedef GenericGeometry::QuadraturePoint< Field, 1 > QuadraturePoint;

      Iterator ( const GaussPoints< Field > &points, const unsigned int index )
        : points_( &points ),
          index_( index )
      {}

      Iterator &operator++ ()
      {
        ++index_;
        return *this;
      }

      QuadraturePoint operator* () const
      {
        return QuadraturePoint( points_->point( index_ ), points_->weight( index_ ) );
      }

    private:
      GaussPoints< Field > *points_;
      unsigned int index_;
    };



    // GaussQuadrature
    // ---------------

    template< class F >
    class GaussQuadrature
      : public Quadrature< 1, F >
    {
      typedef GaussQuadrature< F > This;
      typedef Quadrature< 1, F > Base;

    public:
      typedef typename Base::Field Field;
      static const unsigned int dimension = Base::dimension;

      explicit GaussQuadrature ( unsigned int order )
        : Base( (unsigned int)(0) )
      {
        typedef AlgLib::MultiPrecision< Precision< Field >::value > MPField;
        GaussPoints< MPField > gaussPoints( (order+1) / 2 );
        for( unsigned int i = 0; i < gaussPoints.size(); ++i )
        {
          const Field point = field_cast< Field >( gaussPoints.point( i ) );
          const Field weight = field_cast< Field >( gaussPoints.weight( i ) );
          Base::insert( point, weight );
        }
      }
    };

  }

}

#endif // #ifndef DUNE_GAUSSQUADRATURE_HH
