// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_DUALMORTARBASIS_DUALP1_HH
#define DUNE_LOCALFUNCTIONS_DUALMORTARBASIS_DUALP1_HH

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
   *    Note that if the dual functions are chosen to be dual on the faces,
   *    the integrated product of a Lagrange \f$\lambda_p\f$ and dual
   *    function \f$\theta_q\f$ over faces not containing \f$q\f$ does in
   *    general not vanish.
   *
   * \ingroup DualMortar
   *
   * \tparam D Domain data type
   * \tparam R Range data type
   * \tparam dim Dimension of the simplex
   * \tparam faceDual If set, the basis functions are bi-orthogonal only on faces containing the corresponding vertex.
   */
  template<class D, class R, int dim, bool faceDual=false>
  class DualP1LocalFiniteElement
  {
  public:
    /** \todo Please doc me !
     */
    typedef LocalFiniteElementTraits<DualP1LocalBasis<D,R,dim,faceDual>,DualP1LocalCoefficients<dim>,
        DualP1LocalInterpolation<dim,DualP1LocalBasis<D,R,dim,faceDual> > > Traits;

    /** \todo Please doc me !
     */
    DualP1LocalFiniteElement ()
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

    /** \brief Number of shape functions in this finite element */
    unsigned int size () const
    {
      return basis.size();
    }

    /** \todo Please doc me !
     */
    static constexpr GeometryType type ()
    {
      return GeometryTypes::simplex(dim);
    }

  private:
    DualP1LocalBasis<D,R,dim,faceDual> basis;
    DualP1LocalCoefficients<dim> coefficients;
    DualP1LocalInterpolation<dim,DualP1LocalBasis<D,R,dim,faceDual> > interpolation;
  };



}

#endif
