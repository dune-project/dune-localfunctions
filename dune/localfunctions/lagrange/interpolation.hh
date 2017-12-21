// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LAGRANGEBASIS_INTERPOLATION_HH
#define DUNE_LAGRANGEBASIS_INTERPOLATION_HH

#include <type_traits>
#include <utility>
#include <vector>

#include <dune/common/std/type_traits.hh>
#include <dune/common/typeutilities.hh>

#include <dune/geometry/topologyfactory.hh>

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

    template< class Fn, class Fy >
    auto interpolate ( const Fn &fn, std::vector< Fy > &coefficients, PriorityTag< 1 > ) const
      -> std::enable_if_t< Std::is_invocable< const Fn &, decltype( this->lagrangePoints_.begin()->point() ) >::value >
    {
      unsigned int index = 0;
      for( const auto &lp : lagrangePoints_ )
        field_cast( fn( lp.points() ), coefficients[ index++ ] );
    }

    template< class Fn, class Fy >
    auto interpolate ( const Fn &fn, std::vector< Fy > &coefficients, PriorityTag< 0 > ) const
      -> std::enable_if_t< std::is_same< decltype( fn.evaluate( field_cast< typename Fn::DomainType::field_type >( this->lagrangePoints_.begin()->point() ), std::declval< typename Fn::RangeType & >() ) ), void >::value >
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
    template< class Fn, class Fy >
    void interpolate ( const Fn &fn, std::vector< Fy > &coefficients ) const
    {
      coefficients.resize( lagrangePoints_.size() );
      interpolate( fn, coefficients, PriorityTag< 42 >() );
    }

    template< class Matrix, class Basis >
    void interpolate ( const Basis &basis, Matrix &coefficients ) const
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
  struct LagrangeInterpolationFactoryTraits
  {
    typedef LagrangeCoefficientsFactory<LP,dim,F> LagrangePointSetFactory;
    typedef typename LagrangePointSetFactory::Object LagrangePointSet;

    typedef typename LagrangePointSetFactory::Key Key;
    typedef const LocalLagrangeInterpolation< LP,dim,F > Object;
    typedef LagrangeInterpolationFactory<LP,dim,F> Factory;

    static const unsigned int dimension = dim;
  };

  template< template <class,unsigned int> class LP,
      unsigned int dim, class F >
  struct LagrangeInterpolationFactory :
    public TopologyFactory< LagrangeInterpolationFactoryTraits< LP,dim,F > >
  {
    typedef LagrangeInterpolationFactoryTraits< LP,dim,F > Traits;
    typedef typename Traits::Key Key;
    typedef typename Traits::Object Object;

    template< class Topology >
    static Object *createObject ( const Key &key )
    {
      const typename Traits::LagrangePointSet *lagrangeCoeff
        = Traits::LagrangePointSetFactory::template create< Topology >( key );
      if ( lagrangeCoeff == 0 )
        return 0;
      else
        return new Object( *lagrangeCoeff );
    }
    template< class Topology >
    static bool supports ( const typename Traits::Key &key )
    {
      return true;
    }
    static void release( Object *object)
    {
      Traits::LagrangePointSetFactory::release( object->points() );
      delete object;
    }
  };

}

#endif // #ifndef DUNE_LAGRANGEBASIS_INTERPOLATION_HH
