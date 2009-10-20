// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_L2INTERPOLATION_HH
#define DUNE_L2INTERPOLATION_HH

#include <dune/finiteelements/generic/topologyfactory.hh>
#include <dune/finiteelements/common/matrix.hh>

#include <dune/finiteelements/common/localinterpolation.hh>
#include <dune/finiteelements/quadrature/genericquadrature.hh>

namespace Dune
{
  template< class B, class Q, bool onb >
  struct LocalL2Interpolation;

  template< class B, class Q >
  class LocalL2InterpolationBase
    : public LocalInterpolationInterface< LocalL2InterpolationBase< B, Q > >
  {
    typedef LocalL2InterpolationBase< B, Q > This;

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

  protected:
    LocalL2InterpolationBase ( const Basis &basis, const Quadrature &quadrature )
      : basis_( basis ),
        quadrature_( quadrature )
    {}

    const Basis &basis_;
    const Quadrature &quadrature_;
  };

  template< class B, class Q >
  struct LocalL2Interpolation<B,Q,true>
    : public LocalL2InterpolationBase<B,Q>
  {
    typedef LocalL2InterpolationBase<B,Q> Base;
    template< class BasisFactory, bool onb >
    friend class LocalL2InterpolationFactory;
    using Base::Basis;
    using Base::Quadrature;
  private:
    LocalL2Interpolation ( const typename Base::Basis &basis, const typename Base::Quadrature &quadrature )
      : Base(basis,quadrature)
    {}
  };
  template< class B, class Q >
  struct LocalL2Interpolation<B,Q,false>
    : public LocalL2InterpolationBase<B,Q>
  {
    typedef LocalL2InterpolationBase<B,Q> Base;
    template< class BasisFactory, bool onb >
    friend class LocalL2InterpolationFactory;
    using Base::Basis;
    using Base::Quadrature;
    template< class Function, class DofField >
    void interpolate ( const Function &function, std::vector< DofField > &coefficients ) const
    {
      const unsigned size = Base::basis().size();
      assert( coefficients.size() == size );
      Base::interpolate(function,val_);
      for (unsigned int i=0; i<size; ++i)
      {
        coefficients[i] = 0;
        for (unsigned int j=0; j<size; ++j)
        {
          coefficients[i] += massMatrix_(i,j)*val_[i];
        }
      }
    }
  private:
    LocalL2Interpolation ( const typename Base::Basis &basis, const typename Base::Quadrature &quadrature )
      : Base(basis,quadrature),
        val_(basis.size()),
        massMatrix_()
    {
      typedef typename Base::Quadrature::Iterator Iterator;
      const unsigned size = basis.size();
      massMatrix_.resize( size,size );
      for (unsigned int i=0; i<size; ++i)
        for (unsigned int j=0; j<size; ++j)
          massMatrix_(i,j) = 0; // (i==j)?1:0;
      const Iterator end = quadrature().end();
      for( Iterator it = quadrature().begin(); it != end; ++it )
      {
        basis().evaluate( it->point(), val_ );
        for (unsigned int i=0; i<size; ++i)
          for (unsigned int j=0; j<size; ++j)
            massMatrix_(i,j) += (val_[i]*val_[j])*it->weight(); // (i==j)?1:0;
      }
      massMatrix_.invert();
    }
    typedef typename Base::Basis::Field Field;
    typedef FieldVector< Field, Base::Basis::dimRange > RangeVector;
    typedef LFEMatrix<Field> MassMatrix;
    std::vector<RangeVector> val_;
    MassMatrix massMatrix_;
  };

  template< class BasisFactory, bool onb >
  struct LocalL2InterpolationFactory;
  template< class BasisFactory, bool onb >
  struct LocalL2InterpolationFactoryTraits
  {
    static const unsigned int dimension = BasisFactory::dimension;
    typedef typename BasisFactory::StorageField Field;
    typedef GenericGeometry::Quadrature< dimension, Field > Quadrature;

    typedef typename BasisFactory::Key Key;
    typedef typename BasisFactory::Basis Basis;
    typedef LocalL2Interpolation< Basis, Quadrature, onb > LocalInterpolation;
    typedef const LocalInterpolation Object;
    typedef LocalL2InterpolationFactory<BasisFactory,onb> Factory;
  };

  template< class BasisFactory, bool onb >
  struct LocalL2InterpolationFactory :
    public TopologyFactory< LocalL2InterpolationFactoryTraits<BasisFactory,onb> >
  {
    typedef LocalL2InterpolationFactoryTraits<BasisFactory,onb> Traits;
    static const unsigned int dimension = Traits::dimension;
    typedef typename Traits::Key Key;
    typedef typename Traits::Basis Basis;
    typedef typename Traits::Object Object;;

    typedef typename Traits::Field Field;

    typedef GenericGeometry::Quadrature< dimension, Field > Quadrature;

    template< class Topology >
    static Object *createObject ( const Key &key )
    {
      typedef GenericGeometry::GenericQuadrature< Topology, Field > GenericQuadrature;
      const Basis &basis = *BasisFactory::template create< Topology >( key );
      const Quadrature *quadrature = new GenericQuadrature( 2*basis.order()+1 );
      return new Object( basis, *quadrature );
    }
    static void release ( Object *object )
    {
      const Basis &basis = object->basis();
      const Quadrature &quadrature = object->quadrature();
      BasisFactory::release( &basis );
      delete &quadrature;
      delete object;
    }
  };

}

#endif // #ifndef DUNE_L2INTERPOLATION_HH
