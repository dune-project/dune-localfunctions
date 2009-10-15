// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LABATTOQUADRATURE_HH
#define DUNE_LABATTOQUADRATURE_HH

#include <alglib/gqgenlobatto.h>

#include <dune/common/field.hh>

#include <dune/alglib/multiprecision.hh>
#include <dune/alglib/vector.hh>

#include <dune/finiteelements/quadrature/quadrature.hh>

namespace Dune
{

  namespace GenericGeometry
  {

    // LobattoPoints
    // -----------

    template< class F >
    class LobattoPoints;

    template< unsigned int precision >
    class LobattoPoints< amp::ampf< precision > >
    {
      typedef LobattoPoints< amp::ampf< precision > > This;

    public:
      typedef amp::ampf< precision > Field;
      typedef AlgLib::Vector< Field > Vector;

      struct Iterator;

      explicit LobattoPoints ( unsigned int n )
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



    // LobattoPoints::Iterator
    // ---------------------

    template< unsigned int precision >
    struct LobattoPoints< amp::ampf< precision > >::Iterator
    {
      typedef amp::ampf< precision > Field;
      typedef GenericGeometry::QuadraturePoint< Field, 1 > QuadraturePoint;

      Iterator ( const LobattoPoints< Field > &points, const unsigned int index )
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
      LobattoPoints< Field > *points_;
      unsigned int index_;
    };



    // LobattoQuadrature
    // ---------------

    template< class F >
    class LobattoQuadrature
      : public Quadrature< 1, F >
    {
      typedef LobattoQuadrature< F > This;
      typedef Quadrature< 1, F > Base;

    public:
      typedef typename Base::Field Field;
      static const unsigned int dimension = Base::dimension;

      explicit LobattoQuadrature ( unsigned int order )
        : Base( (unsigned int)(0) )
      {
        typedef amp::ampf< Precision< Field >::value > MPField;
        LobattoPoints< MPField > gaussPoints( (order+1) / 2 );
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
