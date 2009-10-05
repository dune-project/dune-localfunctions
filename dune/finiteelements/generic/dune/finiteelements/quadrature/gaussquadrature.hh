// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GAUSSQUADRATURE_HH
#define DUNE_GAUSSQUADRATURE_HH

#include <alglib/gqgengauss.h>

#include <dune/alglib/multiprecision.hh>
#include <dune/alglib/vector.hh>

#include <dune/grid/quadrature/quadrature.hh>

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

      QuadraturePoint< Field, 1 > operator[] ( const unsigned int i ) const
      {
        return QuadraturePoint< Field, 1 >( point( i ), weight( i ) );
      }

      Field point ( const unsigned int i ) const
      {
        assert( i < points_.size() );
        return points_[ i ];
      }

      Field weight ( const unsigned int i )
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



    // GaussQuadratureRule
    // -------------------

    template< class F >
    class GaussQuadratureRule;

    template< unsigned int precision >
    class GaussQuadratureRule< AlgLib::MultiPrecision< precision > >
      : public QuadratureRule< 1, AlgLib::MultiPrecision< precision > >
    {
      typedef GaussQuadratureRule< AlgLib::MultiPrecision< precision > > This;
      typedef QuadratureRule< 1, AlgLib::MultiPrecision< precision > > Base;

    public:
      typedef typename Base::Field Field;
      static const unsigned int dimension = Base::dimension;

      explicit GaussQuadratureRule ( unsigned int order )
        : Base( GeometryType( GeometryType::cube, 1 ) )
      {
        GaussPoints< Field > gaussPoints( (order+1) / 2 );
        for( unsigned int i = 0; i < gaussPoints.size(); ++i )
          Base::insert( gaussPoints[ i ] );
      }
    };

    template<>
    class GaussQuadratureRule< double >
      : public QuadratureRule< 1, double >
    {
      typedef GaussQuadratureRule< double > This;
      typedef QuadratureRule< 1, double > Base;

    public:
      typedef Base::Field Field;
      static const unsigned int dimension = Base::dimension;

      explicit GaussQuadratureRule ( unsigned int order )
        : Base( GeometryType( GeometryType::cube, 1 ) )
      {
        GaussPoints< AlgLib::MultiPrecision< 128 > > gaussPoints( (order+2) / 2 );
        for( unsigned int i = 0; i < gaussPoints.size(); ++i )
          Base::insert( gaussPoints.point( i ).toDouble(), gaussPoints.weight( i ).toDouble() );
      }
    };
  }

}

#endif // #ifndef DUNE_GAUSSQUADRATURE_HH
