// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_P11DLOCALFINITEELEMENT_HH
#define DUNE_P11DLOCALFINITEELEMENT_HH

#include <dune/common/geometrytype.hh>

#include "common/localfiniteelement.hh"
#include "p11d/p11dlocalbasis.hh"
#include "p11d/p11dlocalcoefficients.hh"
#include "p11d/p11dlocalinterpolation.hh"

namespace Dune
{

  /** \todo Please doc me !
   */
  template<class D, class R>
  class P11DLocalFiniteElement : LocalFiniteElementInterface<
                                     LocalFiniteElementTraits<P11DLocalBasis<D,R>,P11DLocalCoefficients,
                                         P11DLocalInterpolation<P11DLocalBasis<D,R> > >
#ifndef DUNE_VIRTUAL_SHAPEFUNCTIONS
                                     , P11DLocalFiniteElement<D,R>
#endif
                                     >
  {
  public:
    /** \todo Please doc me !
     */
    typedef LocalFiniteElementTraits<P11DLocalBasis<D,R>,P11DLocalCoefficients,
        P11DLocalInterpolation<P11DLocalBasis<D,R> > > Traits;

    /** \todo Please doc me !
     */
    P11DLocalFiniteElement ()
    {
      gt.makeLine();
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
    P11DLocalBasis<D,R> basis;
    P11DLocalCoefficients coefficients;
    P11DLocalInterpolation<P11DLocalBasis<D,R> > interpolation;
    GeometryType gt;
  };

}

#endif
