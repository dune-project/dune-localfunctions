// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_Pk3DLOCALFINITEELEMENT_HH
#define DUNE_Pk3DLOCALFINITEELEMENT_HH

#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include "pk3d/pk3dlocalbasis.hh"
#include "pk3d/pk3dlocalcoefficients.hh"
#include "pk3d/pk3dlocalinterpolation.hh"

namespace Dune
{

  /** \todo Please doc me !
   */
  template<class D, class R, unsigned int k>
  class Pk3DLocalFiniteElement
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

    /** Constructor for variants with permuted vertices.

        \param vertexmap The permutation of the vertices.  This
        can for instance be generated from the global indices of
        the vertices by reducing those to the integers 0...3
     */
    Pk3DLocalFiniteElement (const unsigned int vertexmap[4]) : coefficients(vertexmap)
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

    /** \brief Number of shape functions in this finite element */
    unsigned int size () const
    {
      return basis.size();
    }

    /** \todo Please doc me !
     */
    GeometryType type () const
    {
      return gt;
    }

    Pk3DLocalFiniteElement* clone () const
    {
      return new Pk3DLocalFiniteElement(*this);
    }

  private:
    Pk3DLocalBasis<D,R,k> basis;
    Pk3DLocalCoefficients<k> coefficients;
    Pk3DLocalInterpolation<Pk3DLocalBasis<D,R,k> > interpolation;
    GeometryType gt;
  };

}

#endif
