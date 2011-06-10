// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_Q23DLOCALFINITEELEMENT_HH
#define DUNE_Q23DLOCALFINITEELEMENT_HH

#include <dune/common/geometrytype.hh>

#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include <dune/localfunctions/common/localtoglobaladaptors.hh>
#include "q23d/q23dlocalbasis.hh"
#include "q23d/q23dlocalcoefficients.hh"
#include "q23d/q23dlocalinterpolation.hh"

namespace Dune
{

  /** \todo Please doc me !
   */
  template<class D, class R>
  class Q23DLocalFiniteElement
  {
  public:
    /** \todo Please doc me !
     */
    typedef LocalFiniteElementTraits<Q23DLocalBasis<D,R>,Q23DLocalCoefficients,
        Q23DLocalInterpolation<Q23DLocalBasis<D,R> > > Traits;

    /** \todo Please doc me !
     */
    Q23DLocalFiniteElement ()
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
    Q23DLocalBasis<D,R> basis;
    Q23DLocalCoefficients coefficients;
    Q23DLocalInterpolation<Q23DLocalBasis<D,R> > interpolation;
    GeometryType gt;
  };

  //! Factory for global-valued Q23D elements
  /**
   * \tparam Geometry Type of the geometry.  Used to extract the domain field
   *                  type.
   * \tparam RF       Range field type.
   */
  template<class Geometry, class RF>
  class Q23DFiniteElementFactory :
    public ScalarLocalToGlobalFiniteElementAdaptorFactory<
        Q23DLocalFiniteElement<typename Geometry::ctype, RF>, Geometry
        >
  {
    typedef Q23DLocalFiniteElement<typename Geometry::ctype, RF> LFE;
    typedef ScalarLocalToGlobalFiniteElementAdaptorFactory<LFE, Geometry> Base;

    static const LFE lfe;

  public:
    //! default constructor
    Q23DFiniteElementFactory() : Base(lfe) {}
  };

  template<class Geometry, class RF>
  const typename Q23DFiniteElementFactory<Geometry, RF>::LFE
  Q23DFiniteElementFactory<Geometry, RF>::lfe;
}

#endif
