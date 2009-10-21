// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GENERIC_LOCALFINITEELEMENT_HH
#define DUNE_GENERIC_LOCALFINITEELEMENT_HH

#include <dune/common/typetraits.hh>
#include <dune/common/geometrytype.hh>
#include <dune/grid/genericgeometry/conversion.hh>

#include <dune/finiteelements/common/localfiniteelement.hh>

namespace Dune
{

  // GenericLocalFiniteElementTraits
  // -------------------------------

  template< class BF, class CF, class IF >
  struct GenericLocalFiniteElementTraits
  {
    typedef typename remove_const< typename BF::Object >::type LocalBasisType;
    typedef typename remove_const< typename CF::Object >::type LocalCoefficientsType;
    typedef typename remove_const< typename IF::Object >::type LocalInterpolationType;
  };



  // GenericLocalFiniteElement
  // -------------------------

  template< class BF, class CF, class IF,
      unsigned int dimD, class D, class R >
  class GenericLocalFiniteElement
    : LocalFiniteElementInterface< GenericLocalFiniteElementTraits< BF, CF, IF >,
          GenericLocalFiniteElement< BF, CF, IF, dimD, D, R > >
  {
    typedef GenericLocalFiniteElement< BF, CF, IF, dimD, D, R > This;
    typedef LocalFiniteElementInterface< GenericLocalFiniteElementTraits< BF, CF, IF >, This > Base;

  public:
    typedef typename Base::Traits Traits;

    typedef BF BasisFactory;
    typedef CF CoefficientsFactory;
    typedef IF InterpolationFactory;

    static const int dimDomain = dimD;

    typedef typename BasisFactory::Key Key;

    static_assert( (Conversion< Key, typename CoefficientsFactory::Key >::sameType),
                   "incompatible keys between BasisFactory and CoefficientsFactory" );
    static_assert( (Conversion< Key, typename InterpolationFactory::Key >::sameType),
                   "incompatible keys between BasisFactory and InterpolationFactory" );

    /** \todo Please doc me !
     */
    GenericLocalFiniteElement ( unsigned int topologyId, const Key &key )
      : topologyId_( topologyId ),
        basis_( BasisFactory::create( topologyId, key ) ),
        coefficients_( CoefficientsFactory::create( topologyId, key ) ),
        interpolation_( InterpolationFactory::create( topologyId, key ) )
    {}

    ~GenericLocalFiniteElement()
    {
      BasisFactory::release( basis_ );
      CoefficientsFactory::release( coefficients_ );
      InterpolationFactory::release( interpolation_ );
    }

    /** \todo Please doc me !
     */
    const typename Traits::LocalBasisType &localBasis () const
    {
      if( basis_ == 0 )
        DUNE_THROW( NotImplemented, "LocalBasis not implemented." );
      return *basis_;
    }

    /** \todo Please doc me !
     */
    const typename Traits::LocalCoefficientsType &localCoefficients () const
    {
      if( coefficients_ == 0 )
        DUNE_THROW( NotImplemented, "LocalCoefficients not implemented." );
      return *coefficients_;
    }

    /** \todo Please doc me !
     */
    const typename Traits::LocalInterpolationType &localInterpolation () const
    {
      if( interpolation_ == 0 )
        DUNE_THROW( NotImplemented, "LocalInterpolation not implemented." );
      return *interpolation_;
    }

    /** \todo Please doc me !
     */
    GeometryType type () const
    {
      if( GenericGeometry::hasGeometryType( topologyId_, dimDomain ) )
        return GenericGeometry::geometryType( topologyId_, dimDomain );
      else
        return GeometryType();
    }

    /** \todo Please doc me !
     */
    unsigned int topologyId () const
    {
      return topologyId_;
    }

  private:
    unsigned int topologyId_;
    const typename Traits::LocalBasisType *basis_;
    const typename Traits::LocalCoefficientsType *coefficients_;
    const typename Traits::LocalInterpolationType *interpolation_;
  };

}

#endif // #ifndef DUNE_GENERIC_LOCALFINITEELEMENT_HH
