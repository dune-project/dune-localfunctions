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
    typedef LocalFiniteElementTraits<RefinedP1LocalBasis<D,R,1>,
                                     Impl::LagrangeSimplexLocalCoefficients<1,2>,
                                     Impl::LagrangeSimplexLocalInterpolation<RefinedP1LocalBasis<D,R,1> > > Traits;

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
      return GeometryTypes::line;
    }

  private:
    RefinedP1LocalBasis<D,R,1> basis;
    Impl::LagrangeSimplexLocalCoefficients<1,2> coefficients;
    Impl::LagrangeSimplexLocalInterpolation<RefinedP1LocalBasis<D,R,1> > interpolation;
  };



  /** \todo Please doc me !
   */
  template<class D, class R>
  class RefinedP1LocalFiniteElement<D,R,2>
  {
  public:
    /** \todo Please doc me !
     */
    typedef LocalFiniteElementTraits<RefinedP1LocalBasis<D,R,2>,
                                     Impl::LagrangeSimplexLocalCoefficients<2,2>,
                                     Impl::LagrangeSimplexLocalInterpolation<RefinedP1LocalBasis<D,R,2> > > Traits;

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
      return GeometryTypes::triangle;
    }

  private:
    RefinedP1LocalBasis<D,R,2> basis;
    Impl::LagrangeSimplexLocalCoefficients<2,2> coefficients;
    Impl::LagrangeSimplexLocalInterpolation<RefinedP1LocalBasis<D,R,2> > interpolation;
  };

  /** \todo Please doc me !
   */
  template<class D, class R>
  class RefinedP1LocalFiniteElement<D,R,3>
  {
  public:
    /** \todo Please doc me !
     */
    typedef LocalFiniteElementTraits<RefinedP1LocalBasis<D,R,3>,
                                     Impl::LagrangeSimplexLocalCoefficients<3,2>,
                                     Impl::LagrangeSimplexLocalInterpolation<RefinedP1LocalBasis<D,R,3> > > Traits;

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
      return GeometryTypes::tetrahedron;
    }

  private:
    RefinedP1LocalBasis<D,R,3> basis;
    Impl::LagrangeSimplexLocalCoefficients<3,2> coefficients;
    Impl::LagrangeSimplexLocalInterpolation<RefinedP1LocalBasis<D,R,3> > interpolation;
  };

}

#endif
