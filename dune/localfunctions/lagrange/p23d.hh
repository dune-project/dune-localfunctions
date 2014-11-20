// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_P2_3DLOCALFINITEELEMENT_HH
#define DUNE_P2_3DLOCALFINITEELEMENT_HH

#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include "p23d/p23dlocalbasis.hh"
#include "p23d/p23dlocalcoefficients.hh"
#include "p23d/p23dlocalinterpolation.hh"

namespace Dune
{

  /** \todo Please doc me !
   */
  template<class D, class R>
  class P23DLocalFiniteElement
  {
  public:
    /** \todo Please doc me !
     */
    typedef LocalFiniteElementTraits<P23DLocalBasis<D,R>,
        P23DLocalCoefficients,
        P23DLocalInterpolation<P23DLocalBasis<D,R> > > Traits;

    /** \todo Please doc me !
     */
    P23DLocalFiniteElement ()
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

    P23DLocalFiniteElement* clone () const
    {
      return new P23DLocalFiniteElement(*this);
    }

  private:
    P23DLocalBasis<D,R> basis;
    P23DLocalCoefficients coefficients;
    P23DLocalInterpolation<P23DLocalBasis<D,R> > interpolation;
    GeometryType gt;
  };

}

#endif
