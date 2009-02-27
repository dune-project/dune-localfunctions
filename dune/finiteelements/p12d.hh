// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_P12DLOCALFINITEELEMENT_HH
#define DUNE_P12DLOCALFINITEELEMENT_HH

#include <dune/common/geometrytype.hh>

#include "common/localfiniteelement.hh"
#include "p12d/p12dlocalbasis.hh"
#include "p12d/p12dlocalcoefficients.hh"
#include "p12d/p12dlocalinterpolation.hh"

namespace Dune
{

  /** \todo Please doc me !
   */
  template<class D, class R>
  class P12DLocalFiniteElement : LocalFiniteElementInterface<
                                     LocalFiniteElementTraits<P12DLocalBasis<D,R>,P12DLocalCoefficients,
                                         P12DLocalInterpolation<P12DLocalBasis<D,R> > >,
                                     P12DLocalFiniteElement<D,R> >
  {
  public:
    /** \todo Please doc me !
     */
    typedef LocalFiniteElementTraits<P12DLocalBasis<D,R>,P12DLocalCoefficients,
        P12DLocalInterpolation<P12DLocalBasis<D,R> > > Traits;

    /** \todo Please doc me !
     */
    P12DLocalFiniteElement ()
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
    P12DLocalBasis<D,R> basis;
    P12DLocalCoefficients coefficients;
    P12DLocalInterpolation<P12DLocalBasis<D,R> > interpolation;
    GeometryType gt;
  };

}

#endif
