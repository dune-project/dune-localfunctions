// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_P13DLOCALFINITEELEMENT_HH
#define DUNE_P13DLOCALFINITEELEMENT_HH

#include <dune/common/geometrytype.hh>

#include "common/localfiniteelement.hh"
#include "p13d/p13dlocalbasis.hh"
#include "p13d/p13dlocalcoefficients.hh"
#include "p13d/p13dlocalinterpolation.hh"

namespace Dune
{

  /** \todo Please doc me !
   */
  template<class D, class R>
  class P13DLocalFiniteElement : LocalFiniteElementInterface<
                                     LocalFiniteElementTraits<P13DLocalBasis<D,R>,P13DLocalCoefficients,
                                         P13DLocalInterpolation<P13DLocalBasis<D,R> > >
#ifndef DUNE_VIRTUAL_SHAPEFUNCTIONS
                                     , P13DLocalFiniteElement<D,R>
#endif
                                     >
  {
  public:
    /** \todo Please doc me !
     */
    typedef LocalFiniteElementTraits<P13DLocalBasis<D,R>,P13DLocalCoefficients,
        P13DLocalInterpolation<P13DLocalBasis<D,R> > > Traits;

    /** \todo Please doc me !
     */
    P13DLocalFiniteElement ()
    {
      gt.makeTetrahedron();
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
    P13DLocalBasis<D,R> basis;
    P13DLocalCoefficients coefficients;
    P13DLocalInterpolation<P13DLocalBasis<D,R> > interpolation;
    GeometryType gt;
  };

}

#endif
