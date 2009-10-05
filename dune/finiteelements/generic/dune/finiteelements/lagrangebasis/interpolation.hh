// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LAGRANGEBASIS_INTERPOLATION_HH
#define DUNE_LAGRANGEBASIS_INTERPOLATION_HH

#include <vector>

namespace Dune
{

  template< class LPCreator >
  struct LocalLagrangeInterpolationCreator;

  // LocalLagrangeInterpolation
  // --------------------------

  template< class LPCreator >
  class LocalLagrangeInterpolation
  {
    typedef LocalLagrangeInterpolation< LPCreator > This;

    // template< class LPCreator >
    friend class LocalLagrangeInterpolationCreator< LPCreator >;

  public:
    typedef typename LPCreator::LagrangePoints LagrangePoints;
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



  // LocalLagrangeInterpolationCreator
  // ---------------------------------

  template< class LPCreator >
  struct LocalLagrangeInterpolationCreator
  {
    typedef LPCreator LagrangePointsCreator;
    typedef typename LagrangePointsCreator::Key Key;
    typedef typename LagrangePointsCreator::LagrangePoints LagrangePoints;

    typedef LocalLagrangeInterpolation< LagrangePointsCreator > LocalInterpolation;

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
