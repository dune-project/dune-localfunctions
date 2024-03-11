// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_ENRICHED_SIMPLEXP1BUBBLE_LOCALINTERPOLATION_HH
#define DUNE_LOCALFUNCTIONS_ENRICHED_SIMPLEXP1BUBBLE_LOCALINTERPOLATION_HH

#include <type_traits>
#include <vector>

namespace Dune
{
  /**
   * \ingroup LocalBasisImplementation
   * \brief Interpolation into the SimplexP1BubbleLocalBasis
   *
   * The coefficients `f_i` associated to the vertex DOFs are computed
   * by evaluation of the function `f` to interpolate in the vertex
   * positions `x_i`. The coefficient associated to the bubble function,
   * `f_b`, is obtains by evaluation `f` in the barycenter `x_b` of the
   * element minus the linear contributions in that point:
   *
   * \f[
   *  f_i = f(x_i), i=0,\ldots,dim+1
   *  f_b = f(x_b) - \sum_{i=0}^{dim+1} f_i * \phi_i(x_b),
   * \f]
   *
   * with \f(\phi_i\f) the local basis functions.
   *
   * \nosubgrouping
   **/
  template<class LB>
  class SimplexP1BubbleLocalInterpolation
  {
    static const int dim = LB::dimension;
    static const int numVertices = dim+1;

    using DomainType = typename LB::Traits::DomainType;
    using RangeType = typename LB::Traits::RangeType;

  public:
    /**
     * \brief Local interpolation of the function `f`.
     * \param[in] f  A callable `f:D -> R` with domain `D=DomainType` and range
     *               `R` convertible into the coefficient type `C`.
     * \param[out] out  The interpolation coefficients `{f_i...,f_b}` are stored
     *                  in this output vector.
     **/
    template<class F, class C,
      class R = std::invoke_result_t<F, DomainType>,
      std::enable_if_t<std::is_convertible_v<R, C>, int> = 0>
    static constexpr void interpolate (const F& f, std::vector<C>& out)
    {
      out.resize(numVertices+1);

      // vertices
      DomainType x(0);
      out[0] = f(x);

      for (int i = 0; i < dim; ++i) {
        x = 0;
        x[i] = 1;
        out[i+1] = f(x);
      }

      // element bubble
      x = 1.0/(dim+1);
      R y = f(x);

      // evaluate the other shape functions in x and subtract this value
      std::vector<RangeType> sfValues;
      LB::evaluateFunction(x, sfValues);

      out[numVertices] = y;
      for (int i = 0; i < numVertices; ++i)
        out[numVertices] -= out[i]*sfValues[i];
    }
  };

} // end namespace Dune

#endif // DUNE_LOCALFUNCTIONS_ENRICHED_SIMPLEXP1BUBBLE_LOCALINTERPOLATION_HH
