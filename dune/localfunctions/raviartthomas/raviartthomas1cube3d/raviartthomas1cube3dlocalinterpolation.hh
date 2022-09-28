// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS1_CUBE3D_LOCALINTERPOLATION_HH
#define DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS1_CUBE3D_LOCALINTERPOLATION_HH

#include <vector>

#include <dune/geometry/quadraturerules.hh>
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
  class RT1Cube3DLocalInterpolation
  {

  public:

    /**
     * \brief Make set number s, where 0 <= s < 64
     *
     * \param s Edge orientation indicator
     */
    RT1Cube3DLocalInterpolation (std::bitset<6> s = 0)
    {
      for (size_t i=0; i<6; i++)
        sign_[i] = (s[i]) ? -1.0 : 1.0;

      n_[0] = {-1.0,  0.0,  0.0};
      n_[1] = { 1.0,  0.0,  0.0};
      n_[2] = { 0.0, -1.0,  0.0};
      n_[3] = { 0.0,  1.0,  0.0};
      n_[4] = { 0.0,  0.0, -1.0};
      n_[5] = { 0.0,  0.0,  1.0};
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

      out.resize(36);
      fill(out.begin(), out.end(), 0.0);

      const int qOrder = 3;
      const auto& rule1 = QuadratureRules<Scalar,2>::rule(GeometryTypes::cube(2), qOrder);

      for (auto&& qp : rule1)
      {
        Dune::FieldVector<Scalar,2> qPos = qp.position();
        typename LB::Traits::DomainType localPos;

        localPos = {0.0, qPos[0], qPos[1]};
        auto y = f(localPos);
        out[0] += (y[0]*n_[0][0] + y[1]*n_[0][1] + y[2]*n_[0][2])*qp.weight()*sign_[0];
        out[6] += (y[0]*n_[0][0] + y[1]*n_[0][1] + y[2]*n_[0][2])*(2.0*qPos[0] - 1.0)*qp.weight();
        out[12] += (y[0]*n_[0][0] + y[1]*n_[0][1] + y[2]*n_[0][2])*(2.0*qPos[1] - 1.0)*qp.weight();
        out[18] += (y[0]*n_[0][0] + y[1]*n_[0][1] + y[2]*n_[0][2])*(2.0*qPos[0] - 1.0)*(2.0*qPos[1] - 1.0)*qp.weight();

        localPos = {1.0, qPos[0], qPos[1]};
        y = f(localPos);
        out[1] += (y[0]*n_[1][0] + y[1]*n_[1][1] + y[2]*n_[1][2])*qp.weight()*sign_[1];
        out[7] += (y[0]*n_[1][0] + y[1]*n_[1][1] + y[2]*n_[1][2])*(1.0 - 2.0*qPos[0])*qp.weight();
        out[13] += (y[0]*n_[1][0] + y[1]*n_[1][1] + y[2]*n_[1][2])*(1.0 - 2.0*qPos[1])*qp.weight();
        out[19] += (y[0]*n_[1][0] + y[1]*n_[1][1] + y[2]*n_[1][2])*(1.0 - 2.0*qPos[0])*(2.0*qPos[1] - 1.0)*qp.weight();

        localPos = {qPos[0], 0.0, qPos[1]};
        y = f(localPos);
        out[2] += (y[0]*n_[2][0] + y[1]*n_[2][1] + y[2]*n_[2][2])*qp.weight()*sign_[2];
        out[8] += (y[0]*n_[2][0] + y[1]*n_[2][1] + y[2]*n_[2][2])*(1.0 - 2.0*qPos[0])*qp.weight();
        out[14] += (y[0]*n_[2][0] + y[1]*n_[2][1] + y[2]*n_[2][2])*(2.0*qPos[1] - 1.0)*qp.weight();
        out[20] += (y[0]*n_[2][0] + y[1]*n_[2][1] + y[2]*n_[2][2])*(1.0 - 2.0*qPos[0])*(2.0*qPos[1] - 1.0)*qp.weight();

        localPos = {qPos[0], 1.0, qPos[1]};
        y = f(localPos);
        out[3] += (y[0]*n_[3][0] + y[1]*n_[3][1] + y[2]*n_[3][2])*qp.weight()*sign_[3];
        out[9] += (y[0]*n_[3][0] + y[1]*n_[3][1] + y[2]*n_[3][2])*(2.0*qPos[0] - 1.0)*qp.weight();
        out[15] += (y[0]*n_[3][0] + y[1]*n_[3][1] + y[2]*n_[3][2])*(1.0 - 2.0*qPos[1])*qp.weight();
        out[21] += (y[0]*n_[3][0] + y[1]*n_[3][1] + y[2]*n_[3][2])*(2.0*qPos[0] - 1.0)*(2.0*qPos[1] - 1.0)*qp.weight();

        localPos = {qPos[0], qPos[1], 0.0};
        y = f(localPos);
        out[4] += (y[0]*n_[4][0] + y[1]*n_[4][1] + y[2]*n_[4][2])*qp.weight()*sign_[4];
        out[10] += (y[0]*n_[4][0] + y[1]*n_[4][1] + y[2]*n_[4][2])*(1.0 - 2.0*qPos[0])*qp.weight();
        out[16] += (y[0]*n_[4][0] + y[1]*n_[4][1] + y[2]*n_[4][2])*(1.0 - 2.0*qPos[1])*qp.weight();
        out[22] += (y[0]*n_[4][0] + y[1]*n_[4][1] + y[2]*n_[4][2])*(1.0 - 2.0*qPos[0])*(2.0*qPos[1] - 1.0)*qp.weight();

        localPos = {qPos[0], qPos[1], 1.0};
        y = f(localPos);
        out[5] += (y[0]*n_[5][0] + y[1]*n_[5][1] + y[2]*n_[5][2])*qp.weight()*sign_[5];
        out[11] += (y[0]*n_[5][0] + y[1]*n_[5][1] + y[2]*n_[5][2])*(2.0*qPos[0] - 1.0)*qp.weight();
        out[17] += (y[0]*n_[5][0] + y[1]*n_[5][1] + y[2]*n_[5][2])*(2.0*qPos[1] - 1.0)*qp.weight();
        out[23] += (y[0]*n_[5][0] + y[1]*n_[5][1] + y[2]*n_[5][2])*(2.0*qPos[0] - 1.0)*(2.0*qPos[1] - 1.0)*qp.weight();
      }

      const auto& rule2 = QuadratureRules<Vector,3>::rule(GeometryTypes::cube(3), qOrder);
      for (auto&& qp : rule2)
      {
        FieldVector<double,3> qPos = qp.position();

        auto y = f(qPos);
        out[24] += y[0]*qp.weight();
        out[25] += y[1]*qp.weight();
        out[26] += y[2]*qp.weight();
        out[27] += y[0]*qPos[1]*qp.weight();
        out[28] += y[0]*qPos[2]*qp.weight();
        out[29] += y[1]*qPos[0]*qp.weight();
        out[30] += y[1]*qPos[2]*qp.weight();
        out[31] += y[2]*qPos[0]*qp.weight();
        out[32] += y[2]*qPos[1]*qp.weight();
        out[33] += y[0]*qPos[1]*qPos[2]*qp.weight();
        out[34] += y[1]*qPos[0]*qPos[2]*qp.weight();
        out[35] += y[2]*qPos[0]*qPos[1]*qp.weight();
      }
    }

  private:
    // Facet orientations
    std::array<typename LB::Traits::RangeFieldType, 6> sign_;

    // Facet normals
    std::array<typename LB::Traits::DomainType, 6>     n_;
  };
}
#endif // DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS1_CUBE3D_LOCALINTERPOLATION_HH
