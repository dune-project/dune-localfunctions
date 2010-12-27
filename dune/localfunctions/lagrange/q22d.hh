// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_Q22DLOCALFINITEELEMENT_HH
#define DUNE_Q22DLOCALFINITEELEMENT_HH

#include <dune/common/geometrytype.hh>

#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include <dune/localfunctions/common/localtoglobaladaptors.hh>
#include "q22d/q22dlocalbasis.hh"
#include "q22d/q22dlocalcoefficients.hh"
#include "q22d/q22dlocalinterpolation.hh"

namespace Dune
{

  /** \todo Please doc me !
   */
  template<class D, class R>
  class Q22DLocalFiniteElement
  {
  public:
    /** \todo Please doc me !
     */
    typedef LocalFiniteElementTraits<Q22DLocalBasis<D,R>,Q22DLocalCoefficients,
        Q22DLocalInterpolation<Q22DLocalBasis<D,R> > > Traits;

    /** \todo Please doc me !
     */
    Q22DLocalFiniteElement ()
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
    Q22DLocalBasis<D,R> basis;
    Q22DLocalCoefficients coefficients;
    Q22DLocalInterpolation<Q22DLocalBasis<D,R> > interpolation;
    GeometryType gt;
  };

  //! Factory for global-valued Q22D elements
  /**
   * \tparam Geometry Type of the geometry.  Used to extract the domain field
   *                  type.
   * \tparam RF       Range field type.
   */
  template<class Geometry, class RF>
  class Q22DFiniteElementFactory :
    public ScalarLocalToGlobalFiniteElementAdaptorFactory<
        Q22DLocalFiniteElement<typename Geometry::ctype, RF>, Geometry
        >
  {
    typedef Q22DLocalFiniteElement<typename Geometry::ctype, RF> LFE;
    typedef ScalarLocalToGlobalFiniteElementAdaptorFactory<LFE, Geometry> Base;

    static const LFE lfe;

  public:
    //! default constructor
    Q22DFiniteElementFactory() : Base(lfe) {}
  };

  template<class Geometry, class RF>
  const typename Q22DFiniteElementFactory<Geometry, RF>::LFE
  Q22DFiniteElementFactory<Geometry, RF>::lfe;
}

#endif
