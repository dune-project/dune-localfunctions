// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS2_CUBE2D_LOCALINTERPOLATION_HH
#define DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS2_CUBE2D_LOCALINTERPOLATION_HH

#include <vector>

#include <dune/geometry/quadraturerules.hh>
#include <dune/localfunctions/common/localinterpolation.hh>

namespace Dune
{

  /**
   * \ingroup LocalInterpolationImplementation
   * \brief Second order Raviart-Thomas shape functions on the reference triangle.
   *
   * \tparam LB corresponding LocalBasis giving traits
   *
   * \nosubgrouping
   */
  template<class LB>
  class RT2Cube2DLocalInterpolation
  {

  public:

    /**
     * \brief Make set number s, where 0 <= s < 16
     *
     * \param s Edge orientation indicator
     */
    RT2Cube2DLocalInterpolation (std::bitset<4> s = 0)
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
    template<typename F, typename C>
    void interpolate (const F& ff, std::vector<C>& out) const
    {
      // f gives v*outer normal at a point on the edge!
      typedef typename LB::Traits::RangeFieldType Scalar;
      typedef typename LB::Traits::DomainFieldType Vector;

      auto&& f = Impl::makeFunctionWithCallOperator<typename LB::Traits::DomainType>(ff);

      out.resize(24);
      fill(out.begin(), out.end(), 0.0);

      const int qOrder = 6;
      const auto& rule1 = QuadratureRules<Scalar,1>::rule(GeometryTypes::cube(1), qOrder);

      for (auto&& qp : rule1)
      {
        Scalar qPos = qp.position();
        typename LB::Traits::DomainType localPos;

        localPos = {0.0, qPos};
        auto y = f(localPos);
        out[0] += (y[0]*n_[0][0] + y[1]*n_[0][1])*qp.weight()*sign_[0];
        out[1] += (y[0]*n_[0][0] + y[1]*n_[0][1])*(2.0*qPos - 1.0)*qp.weight();
        out[2] += (y[0]*n_[0][0] + y[1]*n_[0][1])*(6.0*qPos*qPos - 6.0*qPos + 1.0)*qp.weight()*sign_[0];

        localPos = {1.0, qPos};
        y = f(localPos);
        out[3] += (y[0]*n_[1][0] + y[1]*n_[1][1])*qp.weight()*sign_[1];
        out[4] += (y[0]*n_[1][0] + y[1]*n_[1][1])*(1.0 - 2.0*qPos)*qp.weight();
        out[5] += (y[0]*n_[1][0] + y[1]*n_[1][1])*(6.0*qPos*qPos - 6.0*qPos + 1.0)*qp.weight()*sign_[1];

        localPos = {qPos, 0.0};
        y = f(localPos);
        out[6] += (y[0]*n_[2][0] + y[1]*n_[2][1])*qp.weight()*sign_[2];
        out[7] += (y[0]*n_[2][0] + y[1]*n_[2][1])*(1.0 - 2.0*qPos)*qp.weight();
        out[8] += (y[0]*n_[2][0] + y[1]*n_[2][1])*(6.0*qPos*qPos - 6.0*qPos + 1.0)*qp.weight()*sign_[2];

        localPos = {qPos, 1.0};
        y = f(localPos);
        out[9]  += (y[0]*n_[3][0] + y[1]*n_[3][1])*qp.weight()*sign_[3];
        out[10] += (y[0]*n_[3][0] + y[1]*n_[3][1])*(2.0*qPos - 1.0)*qp.weight();
        out[11] += (y[0]*n_[3][0] + y[1]*n_[3][1])*(6.0*qPos*qPos - 6.0*qPos + 1.0)*qp.weight()*sign_[3];
      }

      const auto& rule2 = QuadratureRules<Vector,2>::rule(GeometryTypes::cube(2), qOrder);

      for (auto&& qp : rule2)
      {
        FieldVector<double,2> qPos = qp.position();

        auto y = f(qPos);
        out[12] += y[0]*qp.weight();
        out[13] += y[1]*qp.weight();
        out[14] += y[0]*qPos[0]*qp.weight();
        out[15] += y[1]*qPos[0]*qp.weight();
        out[16] += y[0]*qPos[1]*qp.weight();
        out[17] += y[1]*qPos[1]*qp.weight();
        out[18] += y[0]*qPos[0]*qPos[1]*qp.weight();
        out[19] += y[1]*qPos[0]*qPos[1]*qp.weight();
        out[20] += y[0]*qPos[1]*qPos[1]*qp.weight();
        out[21] += y[1]*qPos[0]*qPos[0]*qp.weight();
        out[22] += y[0]*qPos[0]*qPos[1]*qPos[1]*qp.weight();
        out[23] += y[1]*qPos[0]*qPos[0]*qPos[1]*qp.weight();
      }
    }

  private:
    // Edge orientations
    std::array<typename LB::Traits::RangeFieldType, 4> sign_;

    // Edge normals
    std::array<typename LB::Traits::DomainType, 4>     n_;
  };
}
#endif // DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS2_CUBE2D_LOCALINTERPOLATION_HH
