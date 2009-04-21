// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_Pk3DLOCALFINITEELEMENT_HH
#define DUNE_Pk3DLOCALFINITEELEMENT_HH

#include <dune/common/geometrytype.hh>

#include "common/localfiniteelement.hh"
#include "pk3d/pk3dlocalbasis.hh"
#include "pk3d/pk3dlocalcoefficients.hh"
#include "pk3d/pk3dlocalinterpolation.hh"

namespace Dune
{

  /** \todo Please doc me !
   */
  template<class D, class R, unsigned int k>
  class Pk3DLocalFiniteElement
    :
#if DUNE_VIRTUAL_SHAPEFUNCTIONS
      public
#endif
      LocalFiniteElementInterface<
          LocalFiniteElementTraits<Pk3DLocalBasis<D,R,k>,
              Pk3DLocalCoefficients<k>,
              Pk3DLocalInterpolation<Pk3DLocalBasis<D,R,k> > >
#ifndef DUNE_VIRTUAL_SHAPEFUNCTIONS
          ,Pk3DLocalFiniteElement<D,R,k>
#endif
          >
  {
  public:
    /** \todo Please doc me !
     */
    typedef LocalFiniteElementTraits<Pk3DLocalBasis<D,R,k>,
        Pk3DLocalCoefficients<k>,
        Pk3DLocalInterpolation<Pk3DLocalBasis<D,R,k> > > Traits;

    /** \todo Please doc me !
     */
    Pk3DLocalFiniteElement ()
    {
      gt.makeTetrahedron();
    }

    /** \todo Please doc me !
     */
    Pk3DLocalFiniteElement (unsigned int vertexmap[4]) : coefficients(vertexmap)
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
    Pk3DLocalBasis<D,R,k> basis;
    Pk3DLocalCoefficients<k> coefficients;
    Pk3DLocalInterpolation<Pk3DLocalBasis<D,R,k> > interpolation;
    GeometryType gt;
  };

}

#endif
