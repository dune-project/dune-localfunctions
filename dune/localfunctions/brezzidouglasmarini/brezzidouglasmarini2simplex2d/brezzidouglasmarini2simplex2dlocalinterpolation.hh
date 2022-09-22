// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_BREZZIDOUGLASMARINI2_SIMPLEX2D_LOCALINTERPOLATION_HH
#define DUNE_LOCALFUNCTIONS_BREZZIDOUGLASMARINI2_SIMPLEX2D_LOCALINTERPOLATION_HH

#include <vector>

#include <dune/geometry/quadraturerules.hh>
#include <dune/localfunctions/common/localinterpolation.hh>

namespace Dune
{

  /**
   * \brief First order Brezzi-Douglas-Marini shape functions on triangles.
   *
   * \tparam LB corresponding LocalBasis giving traits
   *
   * \ingroup LocalInterpolationImplementation
   * \nosubgrouping
   */
  template<class LB>
  class BDM2Simplex2DLocalInterpolation
  {

  public:
    //! \brief Standard constructor
    BDM2Simplex2DLocalInterpolation()
    {
      sign0 = sign1 = sign2 = 1.0;
    }

    /**
     * \brief Make set number s, where 0 <= s < 8
     *
     * \param s Edge orientation indicator
     */
    BDM2Simplex2DLocalInterpolation(unsigned int s)
    {
      sign0 = sign1 = sign2 = 1.0;
      if (s & 1)
      {
        sign0 = -1.0;
      }
      if (s & 2)
      {
        sign1 = -1.0;
      }
      if (s & 4)
      {
        sign2 = -1.0;
      }

      m0[0] = 0.5;
      m0[1] = 0.0;
      m1[0] = 0.0;
      m1[1] = 0.5;
      m2[0] = 0.5;
      m2[1] = 0.5;
      n0[0] = 0.0;
      n0[1] = -1.0;
      n1[0] = -1.0;
      n1[1] = 0.0;
      n2[0] = 1.0/sqrt(2.0);
      n2[1] = 1.0/sqrt(2.0);
      c0 =  0.5*n0[0] - 1.0*n0[1];
      c1 = -1.0*n1[0] + 0.5*n1[1];
      c2 =  0.5*n2[0] + 0.5*n2[1];
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
    void interpolate(const F& ff, std::vector<C>& out) const
    {
      // f gives v*outer normal at a point on the edge!
      typedef typename LB::Traits::RangeFieldType Scalar;
      typedef typename LB::Traits::DomainFieldType Vector;

      auto&& f = Impl::makeFunctionWithCallOperator<typename LB::Traits::DomainType>(ff);

      out.resize(12);
      fill(out.begin(), out.end(), 0.0);

      const int qOrder = 4;
      const Dune::QuadratureRule<Scalar,1>& rule = Dune::QuadratureRules<Scalar,1>::rule(Dune::GeometryTypes::simplex(1), qOrder);

      for (typename Dune::QuadratureRule<Scalar,1>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
      {
        Scalar qPos = it->position();

        typename LB::Traits::DomainType localPos;

        localPos[0] = qPos;
        localPos[1] = 0.0;
        auto y = f(localPos);
        out[0] += (y[0]*n0[0] + y[1]*n0[1])*it->weight()*sign0/c0;
        out[1] += (y[0]*n0[0] + y[1]*n0[1])*(1.0 - 2.0*qPos)*it->weight()/c0;
        out[2] += (y[0]*n0[0] + y[1]*n0[1])*(6.0*qPos*qPos - 6.0*qPos + 1.0)*it->weight()*sign0/c0;

        localPos[0] = 0.0;
        localPos[1] = qPos;
        y = f(localPos);
        out[3] += (y[0]*n1[0]+y[1]*n1[1])*it->weight()*sign1/c1;
        out[4] += (y[0]*n1[0]+y[1]*n1[1])*(2.0*qPos-1.0)*it->weight()/c1;
        out[5] += (y[0]*n1[0]+y[1]*n1[1])*(6.0*qPos*qPos - 6.0*qPos + 1.0)*it->weight()*sign1/c1;

        localPos[0] = 1.0 - qPos;
        localPos[1] = qPos;
        y = f(localPos);
        out[6] += (y[0]*n2[0] + y[1]*n2[1])*it->weight()*sign2/c2;
        out[7] += (y[0]*n2[0] + y[1]*n2[1])*(1.0 - 2.0*qPos)*it->weight()/c2;
        out[8] += (y[0]*n2[0] + y[1]*n2[1])*(6.0*qPos*qPos - 6.0*qPos + 1.0)*it->weight()*sign2/c2;
      }

      // a volume part is needed here for dofs: 9 10 11
      const QuadratureRule<Vector,2>& rule2 = QuadratureRules<Vector,2>::rule(GeometryTypes::simplex(2), qOrder);

      for (typename QuadratureRule<Vector,2>::const_iterator it=rule2.begin(); it!=rule2.end(); ++it)
      {
        typename LB::Traits::DomainType localPos = it->position();
        auto y = f(localPos);

        out[9] += y[0]*it->weight();
        out[10] += y[1]*it->weight();
        out[11] += (y[0]*(localPos[0]-2.0*localPos[0]*localPos[1]-localPos[0]*localPos[0])
            +y[1]*(-localPos[1]+2.0*localPos[0]*localPos[1]+localPos[1]*localPos[1]))*it->weight();
      }
    }

  private:
    typename LB::Traits::RangeFieldType sign0, sign1, sign2;
    typename LB::Traits::DomainType m0, m1, m2;
    typename LB::Traits::DomainType n0, n1, n2;
    typename LB::Traits::RangeFieldType c0, c1, c2;
  };
} // end namespace Dune
#endif // DUNE_LOCALFUNCTIONS_BREZZIDOUGLASMARINI2_SIMPLEX2D_LOCALINTERPOLATION_HH
