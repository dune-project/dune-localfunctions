// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DUAL_P1LOCALFINITEELEMENT_HH
#define DUNE_DUAL_P1LOCALFINITEELEMENT_HH

#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include "dualp1/dualp1localbasis.hh"
#include "dualp1/dualp1localcoefficients.hh"
#include "dualp1/dualp1localinterpolation.hh"

namespace Dune
{

  /**
   * \brief The local dual p1 finite element on simplices
   *
   * \ingroup DualMortar
   *
   * \tparam D Domain data type
   * \tparam R Range data type
   * \tparam dim Dimension of the simplex
   */
  template<class D, class R, int dim>
  class DualP1LocalFiniteElement
  {
  public:
    /** \todo Please doc me !
     */
    typedef LocalFiniteElementTraits<DualP1LocalBasis<D,R,dim>,DualP1LocalCoefficients<dim>,
        DualP1LocalInterpolation<dim,DualP1LocalBasis<D,R,dim> > > Traits;

    /** \todo Please doc me !
     */
    DualP1LocalFiniteElement ()
    {
      gt.makeSimplex(dim);
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

    DualP1LocalFiniteElement* clone () const
    {
      return new DualP1LocalFiniteElement(*this);
    }

  private:
    DualP1LocalBasis<D,R,dim> basis;
    DualP1LocalCoefficients<dim> coefficients;
    DualP1LocalInterpolation<dim,DualP1LocalBasis<D,R,dim> > interpolation;
    GeometryType gt;
  };



}

#endif
