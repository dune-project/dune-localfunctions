// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_EDGES02DLOCALFINITEELEMENT_HH
#define DUNE_EDGES02DLOCALFINITEELEMENT_HH

#include <dune/common/geometrytype.hh>

#include "common/localfiniteelementtraits.hh"
#include "edges02d/edges02dlocalbasis.hh"
#include "edges02d/edges02dlocalcoefficients.hh"
#include "edges02d/edges02dlocalinterpolation.hh"

namespace Dune
{

  /** \brief Lowest order 2D edge elements for triangles
   */
  template<class D, class R>
  class EdgeS02DLocalFiniteElement
  {
  public:
    /** \todo Please doc me !
     */
    typedef LocalFiniteElementTraits<
        EdgeS02DLocalBasis<D,R>,
        EdgeS02DLocalCoefficients,
        EdgeS02DLocalInterpolation<EdgeS02DLocalBasis<D,R> >
        > Traits;

    /** \todo Please doc me !
     */
    //! contruct a local finite element instance with default orientations
    EdgeS02DLocalFiniteElement ()
    {
      gt.makeTriangle();
    }

    //! contruct a local finite element instance with the given orientations
    //! \param orientations Bit-map of orientations for each shape function;
    //! bit 0 = 0 means default orientation for the first shape function, bit
    //! 0 = 1 means inverted orientation for the first shape function.
    EdgeS02DLocalFiniteElement (unsigned int s)
      : basis(s), interpolation(s)
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
    EdgeS02DLocalBasis<D,R> basis;
    EdgeS02DLocalCoefficients coefficients;
    EdgeS02DLocalInterpolation<EdgeS02DLocalBasis<D,R> > interpolation;
    GeometryType gt;
  };

}

#endif // DUNE_EDGES02DLOCALFINITEELEMENT_HH
