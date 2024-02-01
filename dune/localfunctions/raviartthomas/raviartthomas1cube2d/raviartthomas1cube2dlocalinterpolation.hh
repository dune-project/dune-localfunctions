// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS1_CUBE2D_LOCALINTERPOLATION_HH
#define DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS1_CUBE2D_LOCALINTERPOLATION_HH

#include <vector>

#include <dune/geometry/quadraturerules.hh>
#include <dune/localfunctions/common/localinterpolation.hh>


namespace Dune
{

  /**
   * \brief First order Raviart-Thomas shape functions on the reference quadrilateral.
   *
   * \tparam LB corresponding LocalBasis giving traits
   *
   * \nosubgrouping
   * \ingroup RaviartThomasImpl
   */
  template<class LB>
  class RT1Cube2DLocalInterpolation
  {

  public:
    /**
     * \brief Make set number s, where 0 <= s < 16
     *
     * \param s Edge orientation indicator
     */
    RT1Cube2DLocalInterpolation (std::bitset<4> s = 0)
    {
      for (size_t i=0; i<4; i++)
        sign_[i] = (s[i]) ? -1.0 : 1.0;

      n_[0] = {-1.0,  0.0};
      n_[1] = { 1.0,  0.0};
      n_[2] = { 0.0, -1.0};
      n_[3] = { 0.0,  1.0};
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
      // f gives v*outer normal at a point on the edge!
      typedef typename LB::Traits::RangeFieldType Scalar;
      typedef typename LB::Traits::DomainFieldType Vector;

      auto&& f = Impl::makeFunctionWithCallOperator<typename LB::Traits::DomainType>(ff);

      out.resize(12);
      fill(out.begin(), out.end(), 0.0);

      const int qOrder = 3;
      const auto& rule1 = QuadratureRules<Scalar,1>::rule(GeometryTypes::cube(1), qOrder);

      for (auto&& qp : rule1)
      {
        Scalar qPos = qp.position();
        typename LB::Traits::DomainType localPos = {0.0, qPos};

        auto y = f(localPos);
        out[0] += (y[0]*n_[0][0] + y[1]*n_[0][1])*qp.weight()*sign_[0];
        out[1] += (y[0]*n_[0][0] + y[1]*n_[0][1])*(2.0*qPos - 1.0)*qp.weight();

        localPos = {1.0, qPos};
        y = f(localPos);
        out[2] += (y[0]*n_[1][0] + y[1]*n_[1][1])*qp.weight()*sign_[1];
        out[3] += (y[0]*n_[1][0] + y[1]*n_[1][1])*(1.0 - 2.0*qPos)*qp.weight();

        localPos = {qPos, 0.0};
        y = f(localPos);
        out[4] += (y[0]*n_[2][0] + y[1]*n_[2][1])*qp.weight()*sign_[2];
        out[5] += (y[0]*n_[2][0] + y[1]*n_[2][1])*(1.0 - 2.0*qPos)*qp.weight();

        localPos = {qPos, 1.0};
        y = f(localPos);
        out[6] += (y[0]*n_[3][0] + y[1]*n_[3][1])*qp.weight()*sign_[3];
        out[7] += (y[0]*n_[3][0] + y[1]*n_[3][1])*(2.0*qPos - 1.0)*qp.weight();
      }

      const auto& rule2 = QuadratureRules<Vector,2>::rule(GeometryTypes::cube(2), qOrder);

      for (auto&& qp : rule2)
      {
        auto qPos = qp.position();

        auto y = f(qPos);
        out[8] += y[0]*qp.weight();
        out[9] += y[1]*qp.weight();
        out[10] += y[0]*qPos[1]*qp.weight();
        out[11] += y[1]*qPos[0]*qp.weight();
      }
    }

  private:
    // Edge orientations
    std::array<typename LB::Traits::RangeFieldType, 4> sign_;

    // Edge normals
    std::array<typename LB::Traits::DomainType, 4>     n_;
  };
}
#endif // DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS1_CUBE2D_LOCALINTERPOLATION_HH
