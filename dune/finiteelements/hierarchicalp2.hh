// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_HIERARCHICAL_P2_LOCALFINITEELEMENT_HH
#define DUNE_HIERARCHICAL_P2_LOCALFINITEELEMENT_HH

#include <dune/common/geometrytype.hh>

#include "common/localfiniteelement.hh"

#include "hierarchicalp2/hierarchicalsimplexp2localbasis.hh"
#include "hierarchicalp2/hierarchicalsimplexp2localinterpolation.hh"

#include "pk2d/pk2dlocalcoefficients.hh"
#include "pk3d/pk3dlocalcoefficients.hh"

namespace Dune
{

  /** \todo Please doc me !
   */
  template<class D, class R, int dim>
  class HierarchicalP2LocalFiniteElement
#ifdef DUNE_VIRTUAL_SHAPEFUNCTIONS
    : public LocalFiniteElementInterface<D,R,dim>
#else
    : LocalFiniteElementInterface<
          LocalFiniteElementTraits<HierarchicalSimplexP2LocalBasis<D,R,dim>,
              typename Dune::SelectType<dim==2, Pk2DLocalCoefficients<2>, Pk3DLocalCoefficients<2> >::Type,
              HierarchicalSimplexP2LocalInterpolation<HierarchicalSimplexP2LocalBasis<D,R,dim> > >,
          HierarchicalP2LocalFiniteElement<D,R,dim> >
#endif
  {

    dune_static_assert(dim==2 || dim==3, "HierarchicalP2LocalFiniteElement only implemented for dim==2, 3.");

  public:
    /** \todo Please doc me !
     */
    typedef LocalFiniteElementTraits<HierarchicalSimplexP2LocalBasis<D,R,dim>,
        typename Dune::SelectType<dim==2, Pk2DLocalCoefficients<2>, Pk3DLocalCoefficients<2> >::Type,
        HierarchicalSimplexP2LocalInterpolation<HierarchicalSimplexP2LocalBasis<D,R,dim> > > Traits;

    /** \todo Please doc me !
     */
    HierarchicalP2LocalFiniteElement ()
    {
      gt.makeTriangle();
    }

    /** \todo Please doc me !
     */
    const typename Traits::LocalBasisType& localBasis () const
    {
      return basis;
    }

    /** \todo Please doc me !
     */
    const typename Traits::LocalCoefficientsType& localCoefficients () const
    {
      return coefficients;
    }

    /** \todo Please doc me !
     */
    const typename Traits::LocalInterpolationType& localInterpolation () const
    {
      return interpolation;
    }

    /** \todo Please doc me !
     */
    GeometryType type () const
    {
      return gt;
    }

  private:
    HierarchicalSimplexP2LocalBasis<D,R,dim> basis;

    /** \todo Stupid, Pk local coefficients can't be parametrized */
    typename Dune::SelectType<dim==2, Pk2DLocalCoefficients<2>, Pk3DLocalCoefficients<2> >::Type coefficients;

    HierarchicalSimplexP2LocalInterpolation<HierarchicalSimplexP2LocalBasis<D,R,dim> > interpolation;
    GeometryType gt;
  };

}

#endif
