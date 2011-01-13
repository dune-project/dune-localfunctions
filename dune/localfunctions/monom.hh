// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_MONOMLOCALFINITEELEMENT_HH
#define DUNE_MONOMLOCALFINITEELEMENT_HH

#include <dune/common/geometrytype.hh>

#include "common/localfiniteelementtraits.hh"
#include "monom/monomlocalbasis.hh"
#include "monom/monomlocalcoefficients.hh"
#include "monom/monomlocalinterpolation.hh"

namespace Dune
{

  /** Monom basis for discontinuous Galerkin

     \tparam D Type used for coordinates
     \tparam R Type used for shape function values
     \tparam d Dimension of the element
     \tparam p Order of the basis
   */
  template<class D, class R, int d, int p>
  class MonomLocalFiniteElement
  {
    enum { static_size = MonomImp::Size<d,p>::val };

  public:
    /** Traits class
     */
    typedef LocalFiniteElementTraits<
        MonomLocalBasis<D,R,d,p>,
        MonomLocalCoefficients<static_size>,
        MonomLocalInterpolation<MonomLocalBasis<D,R,d,p>,static_size>
        > Traits;

    /** \todo Please doc me !
     */
    MonomLocalFiniteElement (GeometryType::BasicType basicType)
      : basis(), interpolation(basicType, basis), gt(basicType,d)
    {}

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
    MonomLocalBasis<D,R,d,p> basis;
    MonomLocalCoefficients<static_size> coefficients;
    MonomLocalInterpolation<MonomLocalBasis<D,R,d,p>,static_size> interpolation;
    GeometryType gt;
  };

}

#endif // DUNE_MONOMLOCALFINITEELEMENT_HH
