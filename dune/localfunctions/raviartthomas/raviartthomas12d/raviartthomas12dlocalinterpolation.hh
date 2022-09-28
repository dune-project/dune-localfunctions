// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS12DLOCALINTERPOLATION_HH
#define DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS12DLOCALINTERPOLATION_HH

#include <vector>

#include <dune/geometry/quadraturerules.hh>
#include <dune/localfunctions/common/localinterpolation.hh>

namespace Dune
{

  /**
   * @ingroup LocalInterpolationImplementation
   * \brief First order Raviart-Thomas shape functions on the reference quadrilateral.
   *
   * \tparam LB corresponding LocalBasis giving traits
   *
   * \nosubgrouping
   */
  template<class LB>
  class RT12DLocalInterpolation
  {

  public:

    /**
     * \brief Make set number s, where 0 <= s < 8
     *
     * \param s Edge orientation indicator
     */
    RT12DLocalInterpolation (std::bitset<3> s = 0)
    {
      using std::sqrt;
      for (size_t i=0; i<3; i++)
        sign_[i] = (s[i]) ? -1.0 : 1.0;

      n_[0] = { 0.0, -1.0};
      n_[1] = {-1.0,  0.0};
      n_[2] = { 1.0/sqrt(2.0), 1.0/sqrt(2.0)};

      c_ = { 0.5*n_[0][0] - 1.0*n_[0][1],
            -1.0*n_[1][0] + 0.5*n_[1][1],
             0.5*n_[2][0] + 0.5*n_[2][1]};
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

      out.resize(8);
      fill(out.begin(), out.end(), 0.0);

      const int qOrder1 = 4;
      const auto& rule1 = Dune::QuadratureRules<Scalar,1>::rule(Dune::GeometryTypes::simplex(1), qOrder1);

      for (auto&& qp : rule1)
      {
        Scalar qPos = qp.position();
        typename LB::Traits::DomainType localPos;

        localPos = {qPos, 0.0};
        auto y = f(localPos);
        out[0] += (y[0]*n_[0][0] + y[1]*n_[0][1])*qp.weight()*sign_[0]/c_[0];
        out[3] += (y[0]*n_[0][0] + y[1]*n_[0][1])*(2.0*qPos - 1.0)*qp.weight()/c_[0];

        localPos = {0.0, qPos};
        y = f(localPos);
        out[1] += (y[0]*n_[1][0] + y[1]*n_[1][1])*qp.weight()*sign_[1]/c_[1];
        out[4] += (y[0]*n_[1][0] + y[1]*n_[1][1])*(1.0 - 2.0*qPos)*qp.weight()/c_[1];

        localPos = {1.0 - qPos, qPos};
        y = f(localPos);
        out[2] += (y[0]*n_[2][0] + y[1]*n_[2][1])*qp.weight()*sign_[2]/c_[2];
        out[5] += (y[0]*n_[2][0] + y[1]*n_[2][1])*(2.0*qPos - 1.0)*qp.weight()/c_[2];
      }

      const int qOrder2 = 8;
      const auto& rule2 = Dune::QuadratureRules<Vector,2>::rule(Dune::GeometryTypes::simplex(2), qOrder2);

      for (auto&& qp : rule2)
      {
        auto qPos = qp.position();

        auto y = f(qPos);
        out[6] += y[0]*qp.weight();
        out[7] += y[1]*qp.weight();
      }
    }

  private:
    // Edge orientations
    std::array<typename LB::Traits::RangeFieldType, 3> sign_;

    // Edge normals
    std::array<typename LB::Traits::DomainType, 3>     n_;

    std::array<typename LB::Traits::RangeFieldType, 3> c_;
  };
}
#endif // DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS12DLOCALINTERPOLATION_HH
