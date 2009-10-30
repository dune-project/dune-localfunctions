// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LABATTOQUADRATURE_HH
#define DUNE_LABATTOQUADRATURE_HH

#if HAVE_ALGLIB
#include <alglib/gqgenlobatto.h>
#endif

#include <dune/finiteelements/generic/math/field.hh>

#include <dune/finiteelements/generic/math/vector.hh>

#include <dune/finiteelements/generic/quadrature/quadrature.hh>

namespace Dune
{

  namespace GenericGeometry
  {

    // LobattoPoints
    // -----------

    template< class F>
    struct LobattoPoints
      : public PointList< F >
    {
      dune_static_assert(sizeof(F)==0,"Lobatto Points only implemented for ALGLib ampf type");
      typedef PointList< F > Base;
      explicit LobattoPoints ( unsigned int n )
        : Base( n )
      {}
    };

#if HAVE_ALGLIB
    template< unsigned int precision >
    class LobattoPoints< amp::ampf< precision > >
      : public PointList< amp::ampf< precision > >
    {
      typedef PointList< amp::ampf< precision > > Base;
      typedef LobattoPoints< amp::ampf< precision > > This;

      using Base::points_;
      using Base::weights_;
    public:
      typedef amp::ampf< precision > Field;

      explicit LobattoPoints ( unsigned int n )
        : Base( n )
      {
        typename Base::Vec alpha( n );
        typename Base::Vec beta( n );
        for( unsigned int i = 0; i < n; ++i )
        {
          alpha[ i ] = 0;
          beta[ i ] = Field( i*i ) / Field( 4*(i*i)-1 );
        }

        bool succ = gqgenlobatto::generategausslobattoquadrature< precision >( alpha, beta, Field( 2 ), Field(-1),Field(1), n, points_, weights_ );

        if (!succ)
        {
          std::cout << "Problem with computing Lobatto points!" << std::endl;
        }
        const Field half = Field( 1 ) / Field( 2 );
        for( unsigned int i = 0; i < n; ++i )
        {
          points_[ i ] = (points_[ i ] + Field( 1 )) * half;
          weights_[ i ] *= half;
        }
      }
    };
#endif

    // LobattoQuadrature
    // ---------------

    template< class F, class CF = F >
    class LobattoQuadrature
      : public Quadrature< 1, F >
    {
      typedef CF ComputeField;
      typedef LobattoQuadrature< F,CF > This;
      typedef Quadrature< 1, F > Base;

    public:
      typedef typename Base::Field Field;
      static const unsigned int dimension = Base::dimension;

      explicit LobattoQuadrature ( unsigned int order )
        : Base( (unsigned int)(0) )
      {
        // typedef amp::ampf< Precision< Field >::value > MPField;
        LobattoPoints< ComputeField > points( (order+1) / 2 );
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
