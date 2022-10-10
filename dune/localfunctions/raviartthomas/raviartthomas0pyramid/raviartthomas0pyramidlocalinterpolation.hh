// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS0_PYRAMID_LOCALINTERPOLATION_HH
#define DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS0_PYRAMID_LOCALINTERPOLATION_HH

#include <vector>

#include <dune/localfunctions/common/localinterpolation.hh>

namespace Dune
{
  /**
   * \ingroup LocalInterpolationImplementation
   * \brief First order Raviart-Thomas shape functions on the reference hexahedron.
   *
   * \tparam LB corresponding LocalBasis giving traits
   *
   * \nosubgrouping
   */
  template<class LB>
  class RT0PyramidLocalInterpolation
  {

  public:

    /**
     * \brief Make set number s, where 0 <= s < 32
     *
     * \param s Face orientation indicator
     */
    RT0PyramidLocalInterpolation (std::bitset<5> s = 0)
    {
      typedef typename LB::Traits::RangeFieldType Scalar;

      for (size_t i=0; i<5; i++)
        sign[i] = (s[i]) ? -1.0 : 1.0;

      Scalar r = 1/std::sqrt(2);

      n[0] = { 0.0,  0.0, -1.0};
      n[1] = {-1.0,  0.0,  0.0};
      n[2] = {   r,  0.0,    r};
      n[3] = { 0.0, -1.0,  0.0};
      n[4] = { 0.0,    r,    r};

      c[0] = 1.0;
      c[1] = 1/2.0;
      c[2] = 1/2.0 * std::sqrt(2);
      c[3] = 1/2.0;
      c[4] = 1/2.0 * std::sqrt(2);

      m[0] = {   0.5,   0.5,   0.0};
      m[1] = {   0.0, 1/3.0, 1/3.0};
      m[2] = { 2/3.0, 1/3.0, 1/3.0};
      m[3] = { 1/3.0,   0.0, 1/3.0};
      m[4] = { 1/3.0, 2/3.0, 1/3.0};
    }

    /**
     * \brief Interpolate a given function with shape functions
     *
     * \tparam F Function type for function which should be interpolated
     * \tparam C Coefficient type
     * \param ff function which should be interpolated
     * \param out return value, vector of coefficients
     */
    template<class F, class C>
    void interpolate (const F& ff, std::vector<C>& out) const
    {
      auto&& f = Impl::makeFunctionWithCallOperator<typename LB::Traits::DomainType>(ff);

      out.resize(5);
      for(int i=0; i<5; i++)
        out[i] = f(m[i]).dot(n[i]) * c[i] * sign[i];
    }

  private:
    // Facet orientations
    std::array<typename LB::Traits::RangeFieldType, 5> sign;
    // Facet area
    std::array<typename LB::Traits::RangeFieldType, 5> c;

    // Facet normals
    std::array<typename LB::Traits::DomainType, 5> n;
    // Facet midpoints
    std::array<typename LB::Traits::DomainType, 5> m;
  };
}
#endif // DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS0_PYRAMID_LOCALINTERPOLATION_HH
