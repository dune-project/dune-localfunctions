// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LAGRANGEBASIS_INTERPOLATION_HH
#define DUNE_LAGRANGEBASIS_INTERPOLATION_HH

#include <vector>
#include <dune/finiteelements/generic/topologyfactory.hh>
#include <dune/finiteelements/common/localinterpolation.hh>

namespace Dune
{

  template< class LPFactory >
  struct LagrangeInterpolationFactory;

  // LocalLagrangeInterpolation
  // --------------------------

  template< class LPFactory >
  class LocalLagrangeInterpolation
    : public LocalInterpolationInterface< LocalLagrangeInterpolation<LPFactory> >
  {
    typedef LocalLagrangeInterpolation< LPFactory > This;

    // template< class LPFactory >
    friend class LagrangeInterpolationFactory< LPFactory >;

  public:
    typedef typename LPFactory::LagrangePoints LagrangePoints;
    typedef typename LagrangePoints::Field Field;

    static const unsigned int dimension = LagrangePoints::dimension;

  private:
    const LagrangePoints &lagrangePoints_;

    explicit LocalLagrangeInterpolation ( const LagrangePoints &lagrangePoints )
      : lagrangePoints_( lagrangePoints )
    {}

  public:
    template< class Function, class Fy >
    void interpolate ( const Function &function, std::vector< Fy > &coefficients ) const
    {
      typedef typename LagrangePoints::iterator Iterator;

      coefficients.resize( lagrangePoints_.size() );

      unsigned int index = 0;
      const Iterator end = lagrangePoints_.end();
      for( Iterator it = lagrangePoints_.begin(); it != end; ++it )
        field_cast(function( it->point() ), coefficients[ index++ ] );
    }

    template< class Matrix, class Basis >
    void interpolate ( const Basis &basis, Matrix &coefficients ) const
    {
      typedef typename LagrangePoints::iterator Iterator;

      coefficients.resize( lagrangePoints_.size(), basis.size( ) );

      unsigned int index = 0;
      const Iterator end = lagrangePoints_.end();
      for( Iterator it = lagrangePoints_.begin(); it != end; ++it )
        basis.template evaluate<0>( it->point(), coefficients.rowPtr( index++ ) );
    }

    const LagrangePoints &lagrangePoints () const
    {
      return lagrangePoints_;
    }
  };



  // LocalLagrangeInterpolationFactory
  // ---------------------------------
  template< class LPFactory >
  struct LagrangeInterpolationFactoryTraits
  {
    typedef typename LPFactory::Key Key;
    typedef const LocalLagrangeInterpolation< LPFactory > Object;
    typedef LagrangeInterpolationFactory<LPFactory> Factory;
  };

  template< class LPFactory >
  struct LagrangeInterpolationFactory :
    public TopologyFactory< LagrangeInterpolationFactoryTraits< LPFactory > >
  {
    typedef LagrangeInterpolationFactoryTraits< LPFactory > Traits;
    typedef typename Traits::Key Key;
    typedef typename Traits::Object Object;
    typedef typename LPFactory::LagrangePoints LagrangePoints;

    template< class Topology >
    static Object *createObject ( const Key &key )
    {
      const LagrangePoints *lagrangePoints
        = LPFactory::template lagrangePoints< Topology >( key );
      if ( lagrangePoints == 0 )
        return 0;
      else
        return new Object( *lagrangePoints );
    }
  };

}

#endif // #ifndef DUNE_LAGRANGEBASIS_INTERPOLATION_HH
