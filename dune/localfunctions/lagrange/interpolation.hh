// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LAGRANGEBASIS_INTERPOLATION_HH
#define DUNE_LAGRANGEBASIS_INTERPOLATION_HH

#include <vector>
#include <dune/geometry/topologyfactory.hh>
#include <dune/localfunctions/lagrange/lagrangecoefficients.hh>

namespace Dune
{

  template< template <class,unsigned int> class LP,
      unsigned int dim, class F >
  struct LagrangeInterpolationFactory;

  // LocalLagrangeInterpolation
  // --------------------------

  template< template <class,unsigned int> class LP,
      unsigned int dim, class F >
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
    const LagrangePointSet *points () const
    {
      return &lagrangePoints_;
    }

  public:
    template< class Function, class Fy >
    void interpolate ( const Function &function, std::vector< Fy > &coefficients ) const
    {
      typedef typename LagrangePointSet::iterator Iterator;

      coefficients.resize( lagrangePoints_.size() );

      unsigned int index = 0;
      const Iterator end = lagrangePoints_.end();
      for( Iterator it = lagrangePoints_.begin(); it != end; ++it )
      {
        typename Function::RangeType val;
        function.evaluate( field_cast<typename Function::DomainType::field_type>(it->point()), val );
        field_cast( val, coefficients[ index++ ] );
      }
    }

    template< class Matrix, class Basis >
    void interpolate ( const Basis &basis, Matrix &coefficients ) const
    {
      typedef typename LagrangePointSet::iterator Iterator;

      coefficients.resize( lagrangePoints_.size(), basis.size( ) );

      unsigned int index = 0;
      const Iterator end = lagrangePoints_.end();
      for( Iterator it = lagrangePoints_.begin(); it != end; ++it )
        basis.template evaluate<0>( it->point(), coefficients.rowPtr( index++ ) );
    }

    const LagrangePointSet &lagrangePoints () const
    {
      return lagrangePoints_;
    }
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
