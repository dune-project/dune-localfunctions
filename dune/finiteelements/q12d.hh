// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_Q12DLOCALFINITEELEMENT_HH
#define DUNE_Q12DLOCALFINITEELEMENT_HH

#include <dune/common/geometrytype.hh>

#include "common/localfiniteelement.hh"
#include "q12d/q12dlocalbasis.hh"
#include "q12d/q12dlocalcoefficients.hh"
#include "q12d/q12dlocalinterpolation.hh"

namespace Dune
{

  /** \todo Please doc me !
   */
  template<class D, class R>
  class Q12DLocalFiniteElement : public LocalFiniteElementInterface<
                                     LocalFiniteElementTraits<Q12DLocalBasis<D,R>,Q12DLocalCoefficients,
                                         Q12DLocalInterpolation<Q12DLocalBasis<D,R> > >
#ifndef DUNE_VIRTUAL_SHAPEFUNCTIONS
                                     ,Q12DLocalFiniteElement<D,R>
#endif
                                     >
  {
  public:
    /** \todo Please doc me !
     */
    typedef LocalFiniteElementTraits<Q12DLocalBasis<D,R>,Q12DLocalCoefficients,
        Q12DLocalInterpolation<Q12DLocalBasis<D,R> > > Traits;

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
    Q12DLocalBasis<D,R> basis;
    Q12DLocalCoefficients coefficients;
    Q12DLocalInterpolation<Q12DLocalBasis<D,R> > interpolation;
    GeometryType gt;
  };

}

#endif
