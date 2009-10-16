// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GAUSSQUADRATURE_HH
#define DUNE_GAUSSQUADRATURE_HH

#if HAVE_ALGLIB
#include <alglib/gqgengauss.h>
#endif

#include <dune/finiteelements/common/field.hh>

#include <dune/grid/common/quadraturerules.hh>

#include <dune/finiteelements/common/vector.hh>

#include <dune/finiteelements/quadrature/quadrature.hh>

namespace Dune
{

  namespace GenericGeometry
  {

    // GaussPoints
    // -----------

    template< class F>
    class GaussPoints
      : public PointList< F >
    {
      typedef PointList< F > Base;
      typedef GaussPoints< F > This;

      using Base::points_;
      using Base::weights_;
    public:
      explicit GaussPoints ( unsigned int n )
        : Base( n )
      {
        const QuadratureRule<double,1>& points = QuadratureRules<double,1>::rule(GeometryType(GeometryType::cube,1), 2*n-1, QuadratureType::Gauss);
        for( unsigned int i = 0; i < n; ++i )
        {
          points_[ i ] = points[i].position()[0];
          weights_[ i ] = points[i].weight();
        }
      }
    };

#if HAVE_ALGLIB
    template< unsigned int precision >
    class GaussPoints< amp::ampf< precision > >
      : public PointList< amp::ampf< precision > >
    {
      typedef PointList< amp::ampf< precision > > Base;
      typedef GaussPoints< amp::ampf< precision > > This;

      using Base::points_;
      using Base::weights_;
    public:
      typedef amp::ampf< precision > Field;

      explicit GaussPoints ( unsigned int n )
        : Base( n )
      {
        typename Base::Vec alpha( n );
        typename Base::Vec beta( n );
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
    };
#endif

    // GaussQuadrature
    // ---------------

    template< class F, class CF = F >
    class GaussQuadrature
      : public Quadrature< 1, F >
    {
      typedef CF ComputeField;
      typedef GaussQuadrature< F,CF > This;
      typedef Quadrature< 1, F > Base;

    public:
      typedef typename Base::Field Field;
      static const unsigned int dimension = Base::dimension;

      explicit GaussQuadrature ( unsigned int order )
        : Base( (unsigned int)(0) )
      {
        // typedef amp::ampf< Precision< Field >::value > MPField;
        GaussPoints< ComputeField > points( (order+1) / 2 );
        for( unsigned int i = 0; i < points.size(); ++i )
        {
          const Field point = field_cast< Field >( points.point( i ) );
          const Field weight = field_cast< Field >( points.weight( i ) );
          Base::insert( point, weight );
        }
      }
    };

  }

}

#endif // #ifndef DUNE_GAUSSQUADRATURE_HH
