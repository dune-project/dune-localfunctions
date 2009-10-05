// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_Q12DLOCALFINITEELEMENT_HH
#define DUNE_Q12DLOCALFINITEELEMENT_HH

#include <dune/common/geometrytype.hh>

#include "common/localfiniteelement.hh"
#include "q1/q1localbasis.hh"
#include "q1/q1localcoefficients.hh"
#include "q1/q1localinterpolation.hh"

namespace Dune
{

  /** \todo Please doc me !
   */
  template<class D, class R>
  class Q12DLocalFiniteElement : public LocalFiniteElementInterface<
                                     LocalFiniteElementTraits<Q1LocalBasis<D,R,2>,Q1LocalCoefficients<2>,
                                         Q1LocalInterpolation<2,Q1LocalBasis<D,R,2> > >
#ifndef DUNE_VIRTUAL_SHAPEFUNCTIONS
                                     ,Q12DLocalFiniteElement<D,R>
#endif
                                     >
  {
  public:
    /** \todo Please doc me !
     */
    typedef LocalFiniteElementTraits<Q1LocalBasis<D,R,2>,Q1LocalCoefficients<2>,
        Q1LocalInterpolation<2,Q1LocalBasis<D,R,2> > > Traits;

    /** \todo Please doc me !
     */
    Q12DLocalFiniteElement ()
    {
      gt.makeQuadrilateral();
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
    Q1LocalBasis<D,R,2> basis;
    Q1LocalCoefficients<2> coefficients;
    Q1LocalInterpolation<2,Q1LocalBasis<D,R,2> > interpolation;
    GeometryType gt;
  };

}

#endif
