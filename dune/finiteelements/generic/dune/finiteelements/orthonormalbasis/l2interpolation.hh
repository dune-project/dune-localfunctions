// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_L2INTERPOLATION_HH
#define DUNE_L2INTERPOLATION_HH

#include <dune/finiteelements/common/localinterpolation.hh>
#include <dune/finiteelements/quadrature/genericquadrature.hh>

namespace Dune
{

  template< class B, class Q >
  class LocalL2Interpolation
    : public LocalInterpolationInterface< LocalL2Interpolation< B, Q > >
  {
    typedef LocalL2Interpolation< B, Q > This;

    template< class BasisCreator >
    friend class LocalL2InterpolationCreator;

  public:
    typedef B Basis;
    typedef Q Quadrature;

    static const unsigned int dimension = Basis::dimension;

    template< class Function, class DofField >
    void interpolate ( const Function &function, std::vector< DofField > &coefficients ) const
    {
      typedef typename Quadrature::Iterator Iterator;
      typedef FieldVector< DofField, Basis::dimRange > RangeVector;

      const unsigned int size = basis().size();
      static std::vector< RangeVector > basisValues( size );

      coefficients.resize( size );
      for( unsigned int i = 0; i < size; ++i )
        coefficients[ i ] = Zero< DofField >();

      const Iterator end = quadrature().end();
      for( Iterator it = quadrature().begin(); it != end; ++it )
      {
        basis().evaluate( it->point(), basisValues );
        RangeVector factor = field_cast< DofField >( function( it->point() ) );
        factor *= field_cast< DofField >( it->weight() );
        for( unsigned int i = 0; i < size; ++i )
          coefficients[ i ] += factor * basisValues[ i ];
      }
    }

    const Basis &basis () const
    {
      return basis_;
    }

    const Quadrature &quadrature () const
    {
      return quadrature_;
    }

  private:
    LocalL2Interpolation ( const Basis &basis, const Quadrature &quadrature )
      : basis_( basis ),
        quadrature_( quadrature )
    {}

    const Basis &basis_;
    const Quadrature &quadrature_;
  };



  template< class BasisCreator >
  struct LocalL2InterpolationCreator
  {
    typedef typename BasisCreator::Key Key;
    typedef typename BasisCreator::Basis Basis;

    typedef typename BasisCreator::StorageField Field;
    static const unsigned int dimension = BasisCreator::dimension;

    typedef GenericGeometry::Quadrature< dimension, Field > Quadrature;

    typedef LocalL2Interpolation< Basis, Quadrature > LocalInterpolation;

    template< class Topology >
    static const LocalInterpolation &localInterpolation ( const Key &key )
    {
      typedef GenericGeometry::GenericQuadrature< Topology, Field > GenericQuadrature;
      const Basis &basis = BasisCreator::template basis< Topology >( key );
      const Quadrature *quadrature = new GenericQuadrature( 2*basis.order()+1 );
      return *(new LocalInterpolation( basis, *quadrature ));
    }

    static void release ( const LocalInterpolation &localInterpolation )
    {
      const Basis &basis = localInterpolation.basis();
      const Quadrature &quadrature = localInterpolation.quadrature();
      delete &localInterpolation;
      delete &quadrature;
      BasisCreator::release( basis );
    }
  };

}

#endif // #ifndef DUNE_L2INTERPOLATION_HH
