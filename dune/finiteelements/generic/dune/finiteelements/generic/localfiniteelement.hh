// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GENERIC_LOCALFINITEELEMENT_HH
#define DUNE_GENERIC_LOCALFINITEELEMENT_HH

#include <dune/common/geometrytype.hh>
#include <dune/grid/genericgeometry/conversion.hh>

#include <dune/finiteelements/common/localfiniteelement.hh>
#include <dune/finiteelements/generic/basisprovider.hh>

namespace Dune
{
  // forward declaration
  template< class BasisC, class CoeffC, class InterpolC,
      unsigned int dimDomain, class D, class R >
  struct GenericLocalFiniteElement
    : LocalFiniteElementInterface<
          LocalFiniteElementTraits< typename BasisC::Basis,
              typename CoeffC::LocalCoefficients,
              typename InterpolC::LocalInterpolation >,
          GenericLocalFiniteElement<BasisC, CoeffC, InterpolC, dimDomain,D,R> >
  {
    typedef GenericLocalFiniteElement<BasisC, CoeffC, InterpolC, dimDomain,D,R> This;
    typedef LocalFiniteElementInterface<
        LocalFiniteElementTraits<
            typename BasisC::Basis,
            typename CoeffC::LocalCoefficients,
            typename InterpolC::LocalInterpolation >, This > Base;
    using Base::Traits;

    typedef typename BasisC::Key Key;

    static_assert( (Conversion<Key,typename CoeffC::Key>::sameType),
                   "incompatible keys between BasisCreator and CoefficientsCreator" );
    static_assert( (Conversion<Key,typename InterpolC::Key>::sameType),
                   "incompatible keys between BasisCreator and InterpolationCreator" );

    typedef FiniteElementProvider<BasisC,CoeffC,InterpolC> FECreator;
    typedef typename FECreator::FiniteElement FiniteElement;

    /** \todo Please doc me !
     */
    GenericLocalFiniteElement ( unsigned int topologyId,
                                const Key &key )
      : topologyId_(topologyId),
        finiteElement_( FECreator::finiteElement(topologyId,key) )
    {}
    ~GenericLocalFiniteElement()
    {
      FECreator::release( finiteElement_ );
    }

    /** \todo Please doc me !
     */
    const typename Traits::LocalBasisType& localBasis () const
    {
      return finiteElement_.basis();
    }

    /** \todo Please doc me !
     */
    const typename Traits::LocalCoefficientsType& localCoefficients () const
    {
      return finiteElement_.coefficients();
    }

    /** \todo Please doc me !
     */
    const typename Traits::LocalInterpolationType& localInterpolation () const
    {
      return finiteElement_.interpolation();
    }

    /** \todo Please doc me !
     */
    GeometryType type () const
    {
      if ( GenericGeometry::hasGeometryType( topologyId_, dimDomain ) )
        return GenericGeometry::geometryType( topologyId_, dimDomain );
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
    const FiniteElement &finiteElement_;
  };

}

#endif
