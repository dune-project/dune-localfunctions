// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_RANNACHER_TUREK2DLOCALFINITEELEMENT_HH
#define DUNE_RANNACHER_TUREK2DLOCALFINITEELEMENT_HH

#include <dune/common/geometrytype.hh>

#include <dune/localfunctions/common/localfiniteelementtraits.hh>

#include "rannacherturek2d/rannacherturek2dlocalbasis.hh"
#include "rannacherturek2d/rannacherturek2dlocalcoefficients.hh"
#include "rannacherturek2d/rannacherturek2dlocalinterpolation.hh"

namespace Dune {

  template<class D, class R>
  class RannacherTurek2DLocalFiniteElement
  {
    RannacherTurek2DLocalBasis<D,R> basis;
    RannacherTurek2DLocalCoefficients coefficients;
    RannacherTurek2DLocalInterpolation<RannacherTurek2DLocalBasis<D,R> > interpolation;
    GeometryType gt;

  public:
    typedef LocalFiniteElementTraits<
        RannacherTurek2DLocalBasis<D,R>,
        RannacherTurek2DLocalCoefficients,
        RannacherTurek2DLocalInterpolation<RannacherTurek2DLocalBasis<D,R> > > Traits;

    RannacherTurek2DLocalFiniteElement () { gt.makeQuadrilateral(); }

    const typename Traits::LocalBasisType& localBasis () const
    {
      return basis;
    }

    const typename Traits::LocalCoefficientsType& localCoefficients () const
    {
      return coefficients;
    }

    const typename Traits::LocalInterpolationType& localInterpolation () const
    {
      return interpolation;
    }

    GeometryType type () const { return gt; }
  };

} // namespace Dune

#endif // DUNE_RANNACHER_TUREK2DLOCALFINITEELEMENT_HH
