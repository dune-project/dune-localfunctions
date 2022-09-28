// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LAGRANGEBASIS_INTERPOLATION_HH
#define DUNE_LAGRANGEBASIS_INTERPOLATION_HH

#include <type_traits>
#include <utility>
#include <vector>

#include <dune/common/typeutilities.hh>

#include <dune/localfunctions/common/localinterpolation.hh>
#include <dune/localfunctions/lagrange/lagrangecoefficients.hh>

namespace Dune
{

  template< template <class,unsigned int> class LP,
      unsigned int dim, class F >
  struct LagrangeInterpolationFactory;

  // LocalLagrangeInterpolation
  // --------------------------

  template< template <class,unsigned int> class LP, unsigned int dim, class F >
  class LocalLagrangeInterpolation
  {
    typedef LocalLagrangeInterpolation< LP,dim,F > This;

  public:
    typedef LP<F,dim> LagrangePointSet;
    typedef typename LagrangePointSet::Field Field;

    static const unsigned int dimension = LagrangePointSet::dimension;

  private:
    friend struct LagrangeInterpolationFactory<LP,dim,F>;
    const LagrangePointSet &lagrangePoints_;

    explicit LocalLagrangeInterpolation ( const LagrangePointSet &lagrangePoints )
      : lagrangePoints_( lagrangePoints )
    {}

    const LagrangePointSet *points () const { return &lagrangePoints_; }

    template< class Fn, class Vector >
    auto interpolate ( const Fn &fn, Vector &coefficients, PriorityTag< 1 > ) const
      -> std::enable_if_t< std::is_invocable_v< const Fn &, decltype( this->lagrangePoints_.begin()->point() ) > >
    {
      unsigned int index = 0;
      for( const auto &lp : lagrangePoints_ )
        field_cast( fn( lp.point() ), coefficients[ index++ ] );
    }
    template< class Fn, class Vector >
    auto interpolate ( const Fn &fn, Vector &coefficients, PriorityTag< 0 > ) const
       -> std::enable_if_t< models<Impl::FunctionWithEvaluate< typename Fn::DomainType, typename Fn::RangeType >, Fn>(), void>
    {
      unsigned int index = 0;
      for( const auto &lp : lagrangePoints_ )
      {
        typename Fn::RangeType val;
        fn.evaluate( field_cast< typename Fn::DomainType::field_type >( lp.point() ), val );
        field_cast( val, coefficients[ index++ ] );
      }
    }

  public:
    template< class Fn, class Vector >
    auto interpolate ( const Fn &fn, Vector &coefficients ) const
    -> std::enable_if_t< std::is_same< decltype(std::declval<Vector>().resize(1) ),void >::value,void>
    {
      coefficients.resize( lagrangePoints_.size() );
      interpolate( fn, coefficients, PriorityTag< 42 >() );
    }

    template< class Basis, class Matrix >
    auto interpolate ( const Basis &basis, Matrix &coefficients ) const
    -> std::enable_if_t< std::is_same<
           decltype(std::declval<Matrix>().rowPtr(0)), typename Matrix::Field* >::value,void>
    {
      coefficients.resize( lagrangePoints_.size(), basis.size( ) );

      unsigned int index = 0;
      for( const auto &lp : lagrangePoints_ )
        basis.template evaluate< 0 >( lp.point(), coefficients.rowPtr( index++ ) );
    }

    const LagrangePointSet &lagrangePoints () const { return lagrangePoints_; }
  };



  // LocalLagrangeInterpolationFactory
  // ---------------------------------
  template< template <class,unsigned int> class LP,
      unsigned int dim, class F >
  struct LagrangeInterpolationFactory
  {
    typedef LagrangeCoefficientsFactory<LP,dim,F> LagrangePointSetFactory;
    typedef typename LagrangePointSetFactory::Object LagrangePointSet;

    typedef typename LagrangePointSetFactory::Key Key;
    typedef const LocalLagrangeInterpolation< LP,dim,F > Object;

    template< GeometryType::Id geometryId >
    static Object *create ( const Key &key )
    {
      const LagrangePointSet *lagrangeCoeff
        = LagrangePointSetFactory::template create< geometryId >( key );
      if ( lagrangeCoeff == 0 )
        return 0;
      else
        return new Object( *lagrangeCoeff );
    }
    template< GeometryType::Id geometryId >
    static bool supports ( const Key &key )
    {
      return true;
    }
    static void release( Object *object)
    {
      LagrangePointSetFactory::release( object->points() );
      delete object;
    }
  };

}

#endif // #ifndef DUNE_LAGRANGEBASIS_INTERPOLATION_HH
