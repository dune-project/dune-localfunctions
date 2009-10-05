// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LAGRANGEBASIS_INTERPOLATION_HH
#define DUNE_LAGRANGEBASIS_INTERPOLATION_HH

#include <dune/finiteelements/lagrangebasis/lagrangepoints.hh>

namespace Dune
{

  // External Forward Declarations
  // -----------------------------

  template< class Topology, class F >
  class MonomialBasis;



  // LocalLagrangeInterpolation
  // --------------------------

  template< class LP >
  class LocalLagrangeInterpolation
  {
    typedef LocalLagrangeInterpolation< LP > This;

    template< class LagrangePointsCreator >
    friend class LocalLagrangeInterpolationCreator;

  public:
    typedef LP LagrangePoints;
    typedef typename LagrangePoints::Field Field;

    static const unsigned int dimension = LagrangePoints::dimension;

  private:
    const LagrangePoints &lagrangePoints_;

    explicit LocalLagrangeInterpolation ( const LagrangePoints &lagrangePoints )
      : lagrangePoints_( lagrangePoints )
    {}

  public:
    template< class Function >
    void interpolate ( const Function &function, std::vector< Field > &coefficients ) const
    {
      typedef typename LagrangePoints::iterator Iterator;

      coefficients.resize( lagrangePoints_.size() );

      unsigned int index = 0;
      const Iterator end = lagrangePoints_.end();
      for( Iterator it = lagrangePoints_.begin(); it != end; ++it )
        coefficients[ index++ ] = function( it->point() );
    }

    template< class Topology, class Matrix >
    void interpolate ( const MonomialBasis< Topology, Field > &basis, Matrix &coefficients ) const
    {
      typedef typename LagrangePoints::iterator Iterator;

      coefficients.resize( lagrangePoints().size(), basis.size( ) );

      unsigned int index = 0;
      const Iterator end = lagrangePoints().end();
      for( Iterator it = lagrangePoints().begin(); it != end; ++it )
        basis.evaluate( it->point(), coefficients.rowPtr( index++ ) );
    }

    template< class Matrix, class Basis >
    void interpolate ( const Basis &basis, Matrix &coefficients ) const
    {
      typedef typename LagrangePoints::iterator Iterator;

      coefficients.resize( lagrangePoints_.size(), basis.size( ) );

      unsigned int index = 0;
      const Iterator end = lagrangePoints_.end();
      for( Iterator it = lagrangePoints_.begin(); it != end; ++it )
        basis.evaluate( it->point(), coefficients.rowPtr( index++ ) );
    }

    const LagrangePoints &lagrangePoints () const
    {
      return lagrangePoints_;
    }
  };



  // LocalLagrangeInterpolationCreator
  // ---------------------------------

  template< class LagrangePointsCreator >
  struct LocalLagrangeInterpolationCreator
  {
    typedef typename LagrangePointsCreator::Key Key;
    typedef typename LagrangePointsCreator::LagrangePoints LagrangePoints;

    typedef LocalLagrangeInterpolation< LagrangePoints > LocalInterpolation;

    template< class Topology >
    static const LocalInterpolation &localInterpolation ( const Key &key )
    {
      const LagrangePoints &lagrangePoints
        = LagrangePointsCreator::template lagrangePoints< Topology >( key );
      return *(new LocalInterpolation( lagrangePoints ));
    }

    static void release ( const LocalInterpolation &localInterpolation )
    {
      const LagrangePoints &lagrangePoints = localInterpolation.lagrangePoints();
      delete &localInterpolation;
      LagrangePointsCreator::release( lagrangePoints );
    }
  };

}

#endif // #ifndef DUNE_LAGRANGEBASIS_INTERPOLATION_HH
