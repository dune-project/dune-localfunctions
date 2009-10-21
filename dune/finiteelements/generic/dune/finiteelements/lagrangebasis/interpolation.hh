// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LAGRANGEBASIS_INTERPOLATION_HH
#define DUNE_LAGRANGEBASIS_INTERPOLATION_HH

#include <vector>
#include <dune/finiteelements/generic/topologyfactory.hh>
#include <dune/finiteelements/common/localinterpolation.hh>
#include <dune/finiteelements/lagrangebasis/lagrangecoefficients.hh>

namespace Dune
{

  template< template <class,unsigned int> class LC,
      unsigned int dim, class F >
  struct LagrangeInterpolationFactory;

  // LocalLagrangeInterpolation
  // --------------------------

  template< template <class,unsigned int> class LC,
      unsigned int dim, class F >
  class LocalLagrangeInterpolation
    : public LocalInterpolationInterface< LocalLagrangeInterpolation<LC,dim,F> >
  {
    typedef LocalLagrangeInterpolation< LC,dim,F > This;

  public:
    typedef LC<F,dim> LagrangeCoefficients;
    typedef typename LagrangeCoefficients::Field Field;

    static const unsigned int dimension = LagrangeCoefficients::dimension;

  private:
    friend class LagrangeInterpolationFactory<LC,dim,F>;
    const LagrangeCoefficients &lagrangePoints_;

    explicit LocalLagrangeInterpolation ( const LagrangeCoefficients &lagrangePoints )
      : lagrangePoints_( lagrangePoints )
    {}
    const LagrangeCoefficients *points () const
    {
      return &lagrangePoints_;
    }

  public:
    template< class Function, class Fy >
    void interpolate ( const Function &function, std::vector< Fy > &coefficients ) const
    {
      typedef typename LagrangeCoefficients::iterator Iterator;

      coefficients.resize( lagrangePoints_.size() );

      unsigned int index = 0;
      const Iterator end = lagrangePoints_.end();
      for( Iterator it = lagrangePoints_.begin(); it != end; ++it )
        field_cast(function( it->point() ), coefficients[ index++ ] );
    }

    template< class Matrix, class Basis >
    void interpolate ( const Basis &basis, Matrix &coefficients ) const
    {
      typedef typename LagrangeCoefficients::iterator Iterator;

      coefficients.resize( lagrangePoints_.size(), basis.size( ) );

      unsigned int index = 0;
      const Iterator end = lagrangePoints_.end();
      for( Iterator it = lagrangePoints_.begin(); it != end; ++it )
        basis.template evaluate<0>( it->point(), coefficients.rowPtr( index++ ) );
    }

    const LagrangeCoefficients &lagrangePoints () const
    {
      return lagrangePoints_;
    }
  };



  // LocalLagrangeInterpolationFactory
  // ---------------------------------
  template< template <class,unsigned int> class LC,
      unsigned int dim, class F >
  struct LagrangeInterpolationFactoryTraits
  {
    typedef Dune::LagrangeCoefficientsFactory<LC,dim,F> LagrangeCoefficientsFactory;
    typedef typename LagrangeCoefficientsFactory::Object LagrangeCoefficients;

    typedef typename LagrangeCoefficientsFactory::Key Key;
    typedef const LocalLagrangeInterpolation< LC,dim,F > Object;
    typedef LagrangeInterpolationFactory<LC,dim,F> Factory;
  };

  template< template <class,unsigned int> class LC,
      unsigned int dim, class F >
  struct LagrangeInterpolationFactory :
    public TopologyFactory< LagrangeInterpolationFactoryTraits< LC,dim,F > >
  {
    typedef LagrangeInterpolationFactoryTraits< LC,dim,F > Traits;
    typedef typename Traits::Key Key;
    typedef typename Traits::Object Object;

    template< class Topology >
    static Object *createObject ( const Key &key )
    {
      const typename Traits::LagrangeCoefficients *lagrangeCoeff
        = Traits::LagrangeCoefficientsFactory::template create< Topology >( key );
      if ( lagrangeCoeff == 0 )
        return 0;
      else
        return new Object( *lagrangeCoeff );
    }
    static void release( Object *object)
    {
      Traits::LagrangeCoefficientsFactory::release( object->points() );
      delete object;
    }
  };

}

#endif // #ifndef DUNE_LAGRANGEBASIS_INTERPOLATION_HH
