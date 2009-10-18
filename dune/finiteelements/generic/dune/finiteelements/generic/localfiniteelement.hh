// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GENERIC_LOCALFINITEELEMENT_HH
#define DUNE_GENERIC_LOCALFINITEELEMENT_HH

#include <dune/common/geometrytype.hh>
#include <dune/grid/genericgeometry/conversion.hh>

#include <dune/finiteelements/common/localfiniteelement.hh>
#include <dune/finiteelements/generic/basisprovider.hh>
#include <dune/finiteelements/generic/polynomialbasis.hh>

namespace Dune
{

  template< class BasisC, class CoeffC, class InterpolC,
      unsigned int dimDomain, class D, class R, class CF >
  struct GenericLocalFiniteElement
    : LocalFiniteElementInterface<
          LocalFiniteElementTraits< GenericLocalBasis<dimDomain,D,R,typename BasisC::Basis > ,
              typename CoeffC::LocalCoefficients,
              typename InterpolC::LocalInterpolation >,
          GenericLocalFiniteElement<BasisC, CoeffC, InterpolC, dimDomain,D,R,CF> >
  {
    typedef FiniteElementProvider<BasisC,CoeffC,InterpolC> FECreator;
    typedef typename FECreator::FiniteElement FiniteElement;

    /** \todo Please doc me !
     */
    typedef LocalFiniteElementTraits< GenericLocalBasis<dimDomain,D,R,typename FECreator::Basis > ,
        typename FECreator::LocalCoefficients,
        typename FECreator::LocalInterpolation > Traits;

    /** \todo Please doc me !
     */
    GenericLocalFiniteElement ( unsigned int topologyId,
                                unsigned int order )
      : topologyId_(topologyId),
        order_(order),
        finiteElement_( FECreator::finiteElement(topologyId,order) ),
        localBasis_(finiteElement_.basis())
    {}
    ~GenericLocalFiniteElement()
    {
      FECreator::release( finiteElement_ );
    }

    /** \todo Please doc me !
     */
    const typename Traits::LocalBasisType& localBasis () const
    {
      return localBasis_;
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
    unsigned int order_;
    const FiniteElement &finiteElement_;
    typename Traits::LocalBasisType localBasis_;
  };

}

#endif
