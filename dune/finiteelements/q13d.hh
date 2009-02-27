// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_Q13DLOCALFINITEELEMENT_HH
#define DUNE_Q13DLOCALFINITEELEMENT_HH

#include <dune/common/geometrytype.hh>

#include "common/localfiniteelement.hh"
#include "q13d/q13dlocalbasis.hh"
#include "q13d/q13dlocalcoefficients.hh"
#include "q13d/q13dlocalinterpolation.hh"

namespace Dune
{

  /** \todo Please doc me !
   */
  template<class D, class R>
  class Q13DLocalFiniteElement : LocalFiniteElementInterface<
                                     LocalFiniteElementTraits<Q13DLocalBasis<D,R>,Q13DLocalCoefficients,
                                         Q13DLocalInterpolation<Q13DLocalBasis<D,R> > >,
                                     Q13DLocalFiniteElement<D,R> >
  {
  public:
    /** \todo Please doc me !
     */
    typedef LocalFiniteElementTraits<Q13DLocalBasis<D,R>,Q13DLocalCoefficients,
        Q13DLocalInterpolation<Q13DLocalBasis<D,R> > > Traits;

    /** \todo Please doc me !
     */
    Q13DLocalFiniteElement ()
    {
      gt.makeHexahedron();
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
    Q13DLocalBasis<D,R> basis;
    Q13DLocalCoefficients coefficients;
    Q13DLocalInterpolation<Q13DLocalBasis<D,R> > interpolation;
    GeometryType gt;
  };

}

#endif
