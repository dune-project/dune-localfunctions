// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_EDGES12DLOCALFINITEELEMENT_HH
#define DUNE_EDGES12DLOCALFINITEELEMENT_HH

#include <dune/common/geometrytype.hh>

#include "common/localfiniteelement.hh"
#include "edges12d/edges12dlocalbasis.hh"
#include "edges12d/edges12dlocalcoefficients.hh"
#include "edges12d/edges12dlocalinterpolation.hh"

namespace Dune
{

  /** \todo Please doc me !
   */
  template<class D, class R>
  class EdgeS12DLocalFiniteElement
    : public LocalFiniteElementInterface<
          LocalFiniteElementTraits<
              EdgeS12DLocalBasis<D,R>,
              EdgeS12DLocalCoefficients,
              EdgeS12DLocalInterpolation<EdgeS12DLocalBasis<D,R> >
              >
#ifndef DUNE_VIRTUAL_SHAPEFUNCTIONS
          , EdgeS12DLocalFiniteElement<D,R>
#endif
          >
  {
  public:
    /** \todo Please doc me !
     */
    typedef LocalFiniteElementTraits<
        EdgeS12DLocalBasis<D,R>,
        EdgeS12DLocalCoefficients,
        EdgeS12DLocalInterpolation<EdgeS12DLocalBasis<D,R> >
        > Traits;

    /** \todo Please doc me !
     */
    EdgeS12DLocalFiniteElement ()
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
    EdgeS12DLocalBasis<D,R> basis;
    EdgeS12DLocalCoefficients coefficients;
    EdgeS12DLocalInterpolation<EdgeS12DLocalBasis<D,R> > interpolation;
    GeometryType gt;
  };

}

#endif // DUNE_EDGES12DLOCALFINITEELEMENT_HH
