// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_P1LOCALFINITEELEMENT_HH
#define DUNE_P1LOCALFINITEELEMENT_HH

#include <dune/common/geometrytype.hh>

#include "common/localfiniteelement.hh"
#include "p1/p1localbasis.hh"
#include "p1/p1localcoefficients.hh"
#include "p1/p1localinterpolation.hh"

namespace Dune
{

  /** \todo The local p1 finite element on simplices
      \tparam D Domain data type
      \tparam R Range data type
      \tparam dim Dimension of the simplex
   */
  template<class D, class R, int dim>
  class P1LocalFiniteElement
#ifdef DUNE_VIRTUAL_SHAPEFUNCTIONS
    : public LocalFiniteElementInterface<D,R,dim>
#else
    : public LocalFiniteElementInterface<LocalFiniteElementTraits<P1LocalBasis<D,R,dim>,
              P1LocalCoefficients<dim>,
              P1LocalInterpolation<dim,P1LocalBasis<D,R,dim> > >,
          P1LocalFiniteElement<D,R,dim> >
#endif
  {
  public:
    /** \todo Please doc me !
     */
    typedef LocalFiniteElementTraits<P1LocalBasis<D,R,dim>,P1LocalCoefficients<dim>,
        P1LocalInterpolation<dim,P1LocalBasis<D,R,dim> > > Traits;

    /** \todo Please doc me !
     */
    P1LocalFiniteElement ()
    {
      gt.makeSimplex(dim);
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

    P1LocalFiniteElement* clone () const
    {
      return new P1LocalFiniteElement(*this);
    }

  private:
    P1LocalBasis<D,R,dim> basis;
    P1LocalCoefficients<dim> coefficients;
    P1LocalInterpolation<dim,P1LocalBasis<D,R,dim> > interpolation;
    GeometryType gt;
  };



}

#endif
