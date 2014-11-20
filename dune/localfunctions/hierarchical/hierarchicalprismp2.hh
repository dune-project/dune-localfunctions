// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_HIERARCHICAL_PRISM_P2_LOCALFINITEELEMENT_HH
#define DUNE_HIERARCHICAL_PRISM_P2_LOCALFINITEELEMENT_HH

#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include <dune/localfunctions/lagrange/prismp2/prismp2localcoefficients.hh>

#include "hierarchicalprismp2/hierarchicalprismp2localbasis.hh"
#include "hierarchicalprismp2/hierarchicalprismp2localinterpolation.hh"


namespace Dune
{

  /** \todo Please doc me !
   */
  template<class D, class R>
  class HierarchicalPrismP2LocalFiniteElement
  {


  public:
    /** \todo Please doc me !
     */
    typedef LocalFiniteElementTraits<HierarchicalPrismP2LocalBasis<D,R>,
        PrismP2LocalCoefficients,
        HierarchicalPrismP2LocalInterpolation<HierarchicalPrismP2LocalBasis<D,R> > > Traits;

    /** \todo Please doc me !
     */
    HierarchicalPrismP2LocalFiniteElement ()
    {
      gt.makePrism();
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

    /** \brief Number of shape functions in this finite element */
    unsigned int size () const
    {
      return basis.size();
    }

    /** \todo Please doc me !
     */
    GeometryType type () const
    {
      return gt;
    }

    HierarchicalPrismP2LocalFiniteElement* clone () const
    {
      return new HierarchicalPrismP2LocalFiniteElement(*this);
    }

  private:
    HierarchicalPrismP2LocalBasis<D,R> basis;

    /** \todo Stupid, Pk local coefficients can't be parametrized */
    PrismP2LocalCoefficients coefficients;

    HierarchicalPrismP2LocalInterpolation<HierarchicalPrismP2LocalBasis<D,R> > interpolation;
    GeometryType gt;
  };

}

#endif
