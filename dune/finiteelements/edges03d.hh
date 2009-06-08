// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_EDGES03DLOCALFINITEELEMENT_HH
#define DUNE_EDGES03DLOCALFINITEELEMENT_HH

#include <dune/common/geometrytype.hh>

#include "common/localfiniteelement.hh"
#include "edges03d/edges03dlocalbasis.hh"
#include "edges03d/edges03dlocalcoefficients.hh"
#include "edges03d/edges03dlocalinterpolation.hh"

namespace Dune
{

  /** \brief Lowest order 3D edge elements for tetrahedrons
   */
  template<class D, class R>
  class EdgeS03DLocalFiniteElement
    : public LocalFiniteElementInterface<
          LocalFiniteElementTraits<
              EdgeS03DLocalBasis<D,R>,
              EdgeS03DLocalCoefficients,
              EdgeS03DLocalInterpolation<EdgeS03DLocalBasis<D,R> >
              >
#ifndef DUNE_VIRTUAL_SHAPEFUNCTIONS
          , EdgeS03DLocalFiniteElement<D,R>
#endif
          >
  {
  public:
    /** \todo Please doc me !
     */
    typedef LocalFiniteElementTraits<
        EdgeS03DLocalBasis<D,R>,
        EdgeS03DLocalCoefficients,
        EdgeS03DLocalInterpolation<EdgeS03DLocalBasis<D,R> >
        > Traits;

    /** \todo Please doc me !
     */
    //! contruct a local finite element instance with default orientations
    EdgeS03DLocalFiniteElement ()
    {
      gt.makeTetrahedron();
    }

    //! contruct a local finite element instance with the given orientations
    //! \param orientations Bit-map of orientations for each shape function;
    //! bit 0 = 0 means default orientation for the first shape function, bit
    //! 0 = 1 means inverted orientation for the first shape function.
    EdgeS03DLocalFiniteElement (unsigned int s)
      : basis(s), interpolation(s)
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
    EdgeS03DLocalBasis<D,R> basis;
    EdgeS03DLocalCoefficients coefficients;
    EdgeS03DLocalInterpolation<EdgeS03DLocalBasis<D,R> > interpolation;
    GeometryType gt;
  };

}

#endif // DUNE_EDGES03DLOCALFINITEELEMENT_HH
