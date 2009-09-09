// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_HIERARCHICAL_P2_WITH_ELEMENTBUBBLE_LOCALFINITEELEMENT_HH
#define DUNE_HIERARCHICAL_P2_WITH_ELEMENTBUBBLE_LOCALFINITEELEMENT_HH

#include <dune/common/geometrytype.hh>

#include <dune/finiteelements/common/localfiniteelement.hh>
#include <dune/finiteelements/hierarchicalp2withelementbubble/hierarchicalsimplexp2withelementbubble.hh>


namespace Dune
{

  /** \todo Please doc me !
   */
  template<class D, class R, int dim>
  class HierarchicalP2WithElementBubbleLocalFiniteElement
    : LocalFiniteElementInterface<
          LocalFiniteElementTraits<HierarchicalSimplexP2WithElementBubbleLocalBasis<D,R,dim>,
              HierarchicalSimplexP2WithElementBubbleLocalCoefficients<dim>,
              HierarchicalSimplexP2WithElementBubbleLocalInterpolation<HierarchicalSimplexP2WithElementBubbleLocalBasis<D,R,dim> > >
#ifndef DUNE_VIRTUAL_SHAPEFUNCTIONS
          , HierarchicalP2WithElementBubbleLocalFiniteElement<D,R,dim>
#endif
          >
  {

    dune_static_assert(dim==2, "HierarchicalP2WithElementBubbleLocalFiniteElement only implemented for dim==2.");

  public:
    /** \todo Please doc me !
     */
    typedef LocalFiniteElementTraits<HierarchicalSimplexP2WithElementBubbleLocalBasis<D,R,dim>,
        HierarchicalSimplexP2WithElementBubbleLocalCoefficients<dim>,
        HierarchicalSimplexP2WithElementBubbleLocalInterpolation<HierarchicalSimplexP2WithElementBubbleLocalBasis<D,R,dim> > > Traits;

    /** \todo Please doc me !
     */
    HierarchicalP2WithElementBubbleLocalFiniteElement ()
    {
      gt_.makeTriangle();
    }

    /** \todo Please doc me !
     */
    const typename Traits::LocalBasisType& localBasis () const
    {
      return basis_;
    }

    /** \todo Please doc me !
     */
    const typename Traits::LocalCoefficientsType& localCoefficients () const
    {
      return coefficients_;
    }

    /** \todo Please doc me !
     */
    const typename Traits::LocalInterpolationType& localInterpolation () const
    {
      return interpolation_;
    }

    /** \todo Please doc me !
     */
    GeometryType type () const
    {
      return gt_;
    }

  private:
    HierarchicalSimplexP2WithElementBubbleLocalBasis<D,R,dim> basis_;

    HierarchicalSimplexP2WithElementBubbleLocalCoefficients<dim> coefficients_;

    HierarchicalSimplexP2WithElementBubbleLocalInterpolation<HierarchicalSimplexP2WithElementBubbleLocalBasis<D,R,dim> > interpolation_;

    GeometryType gt_;
  };

}

#endif
