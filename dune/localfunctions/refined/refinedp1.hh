// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_REFINED_P1_LOCALFINITEELEMENT_HH
#define DUNE_REFINED_P1_LOCALFINITEELEMENT_HH

#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include <dune/localfunctions/lagrange/p0.hh>

#include <dune/localfunctions/lagrange/lagrangesimplex.hh>
#include <dune/localfunctions/refined/refinedp1/refinedp1localbasis.hh>

namespace Dune
{

  /** \todo Please doc me !
   */
  template<class D, class R, int dim>
  class RefinedP1LocalFiniteElement
  {
  public:
    /** \todo Please doc me !
     */
    typedef LocalFiniteElementTraits<RefinedP1LocalBasis<D,R,dim>,
                                     Impl::LagrangeSimplexLocalCoefficients<dim,2>,
                                     Impl::LagrangeSimplexLocalInterpolation<Impl::LagrangeSimplexLocalBasis<D,R,dim,2> > > Traits;

    /** \todo Please doc me !
     */
    RefinedP1LocalFiniteElement ()
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
      return GeometryTypes::simplex(dim);
    }

  private:
    RefinedP1LocalBasis<D,R,dim> basis;
    Impl::LagrangeSimplexLocalCoefficients<dim,2> coefficients;
    // Yes, the template argument here really is LagrangeSimplexLocalBasis, even though this is not
    // the local basis of the refined locale finite element:  The reason is that LagrangeSimplexLocalInterpolation
    // uses this argument to determine the polynomial order, and RefinedP1LocalBasis returns order 1
    // whereas order 2 is needed here.
    Impl::LagrangeSimplexLocalInterpolation<Impl::LagrangeSimplexLocalBasis<D,R,dim,2> > interpolation;
  };

}

#endif
