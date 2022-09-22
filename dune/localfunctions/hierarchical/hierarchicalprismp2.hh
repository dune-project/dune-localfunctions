// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_HIERARCHICAL_PRISM_P2_LOCALFINITEELEMENT_HH
#define DUNE_HIERARCHICAL_PRISM_P2_LOCALFINITEELEMENT_HH

#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include <dune/localfunctions/lagrange/lagrangeprism.hh>

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
        Impl::LagrangePrismLocalCoefficients<2>,
        HierarchicalPrismP2LocalInterpolation<HierarchicalPrismP2LocalBasis<D,R> > > Traits;

    /** \todo Please doc me !
     */
    HierarchicalPrismP2LocalFiniteElement ()
    {}

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
    static constexpr GeometryType type ()
    {
      return GeometryTypes::prism;
    }

  private:
    HierarchicalPrismP2LocalBasis<D,R> basis;

    /** \todo Stupid, Pk local coefficients can't be parametrized */
    Impl::LagrangePrismLocalCoefficients<2> coefficients;

    HierarchicalPrismP2LocalInterpolation<HierarchicalPrismP2LocalBasis<D,R> > interpolation;
  };

}

#endif
