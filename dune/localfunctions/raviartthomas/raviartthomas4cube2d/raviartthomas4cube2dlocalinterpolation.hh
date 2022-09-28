// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS4_CUBE2D_LOCALINTERPOLATION_HH
#define DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS4_CUBE2D_LOCALINTERPOLATION_HH

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
  class RT4Cube2DLocalInterpolation
  {

  public:

    /**
     * \brief Make set number s, where 0 <= s < 16
     *
     * \param s Edge orientation indicator
     */
    RT4Cube2DLocalInterpolation (unsigned int s = 0)
    {
      sign0 = sign1 = sign2 = sign3 = 1.0;
      if (s & 1)
      {
        sign0 *= -1.0;
      }
      if (s & 2)
      {
        sign1 *= -1.0;
      }
      if (s & 4)
      {
        sign2 *= -1.0;
      }
      if (s & 8)
      {
        sign3 *= -1.0;
      }

      n0[0] = -1.0;
      n0[1] =  0.0;
      n1[0] =  1.0;
      n1[1] =  0.0;
      n2[0] =  0.0;
      n2[1] = -1.0;
      n3[0] =  0.0;
      n3[1] =  1.0;
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

      out.resize(60);
      fill(out.begin(), out.end(), 0.0);

      const int qOrder = 12;
      const QuadratureRule<Scalar,1>& rule = QuadratureRules<Scalar,1>::rule(GeometryTypes::cube(1), qOrder);

      for (typename QuadratureRule<Scalar,1>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
      {
        Scalar qPos = it->position();
        typename LB::Traits::DomainType localPos;

        localPos[0] = 0.0;
        localPos[1] = qPos;
        auto y = f(localPos);
        out[0] += (y[0]*n0[0] + y[1]*n0[1])*it->weight()*sign0;
        out[1] += (y[0]*n0[0] + y[1]*n0[1])*(2.0*qPos - 1.0)*it->weight();
        out[2] += (y[0]*n0[0] + y[1]*n0[1])*(6.0*qPos*qPos - 6.0*qPos + 1.0)*it->weight()*sign0;
        out[3] += (y[0]*n0[0] + y[1]*n0[1])*(20.0*qPos*qPos*qPos - 30.0*qPos*qPos + 12.0*qPos - 1.0)*it->weight();
        out[4] += (y[0]*n0[0] + y[1]*n0[1])*(1.0-20.0*qPos+90.0*pow(qPos,2)-140.0*pow(qPos,3)+70.0*pow(qPos,4))*it->weight()*sign0;

        localPos[0] = 1.0;
        localPos[1] = qPos;
        y = f(localPos);
        out[5] += (y[0]*n1[0] + y[1]*n1[1])*it->weight()*sign1;
        out[6] += (y[0]*n1[0] + y[1]*n1[1])*(1.0 - 2.0*qPos)*it->weight();
        out[7] += (y[0]*n1[0] + y[1]*n1[1])*(6.0*qPos*qPos - 6.0*qPos + 1.0)*it->weight()*sign1;
        out[8] += (y[0]*n1[0] + y[1]*n1[1])*(-20.0*qPos*qPos*qPos + 30.0*qPos*qPos - 12.0*qPos + 1.0)*it->weight();
        out[9] += (y[0]*n1[0] + y[1]*n1[1])*(1.0-20.0*qPos+90.0*pow(qPos,2)-140.0*pow(qPos,3)+70.0*pow(qPos,4))*it->weight()*sign1;

        localPos[0] = qPos;
        localPos[1] = 0.0;
        y = f(localPos);
        out[10] += (y[0]*n2[0] + y[1]*n2[1])*it->weight()*sign2;
        out[11] += (y[0]*n2[0] + y[1]*n2[1])*(1.0 - 2.0*qPos)*it->weight();
        out[12] += (y[0]*n2[0] + y[1]*n2[1])*(6.0*qPos*qPos - 6.0*qPos + 1.0)*it->weight()*sign2;
        out[13] += (y[0]*n2[0] + y[1]*n2[1])*(-20.0*qPos*qPos*qPos + 30.0*qPos*qPos - 12.0*qPos + 1.0)*it->weight();
        out[14] += (y[0]*n2[0] + y[1]*n2[1])*(1.0-20.0*qPos+90.0*pow(qPos,2)-140.0*pow(qPos,3)+70.0*pow(qPos,4))*it->weight()*sign2;

        localPos[0] = qPos;
        localPos[1] = 1.0;
        y = f(localPos);
        out[15]  += (y[0]*n3[0] + y[1]*n3[1])*it->weight()*sign3;
        out[16] += (y[0]*n3[0] + y[1]*n3[1])*(2.0*qPos - 1.0)*it->weight();
        out[17] += (y[0]*n3[0] + y[1]*n3[1])*(6.0*qPos*qPos - 6.0*qPos + 1.0)*it->weight()*sign3;
        out[18] += (y[0]*n3[0] + y[1]*n3[1])*(20.0*qPos*qPos*qPos - 30.0*qPos*qPos + 12.0*qPos - 1.0)*it->weight();
        out[19] += (y[0]*n3[0] + y[1]*n3[1])*(1.0-20.0*qPos+90.0*pow(qPos,2)-140.0*pow(qPos,3)+70.0*pow(qPos,4))*it->weight()*sign3;
      }

      const QuadratureRule<Vector,2>& rule2 = QuadratureRules<Vector,2>::rule(GeometryTypes::cube(2), qOrder);

      for (typename QuadratureRule<Vector,2>::const_iterator it = rule2.begin();
           it != rule2.end(); ++it)
      {
        FieldVector<double,2> qPos = it->position();

        auto y = f(qPos);
        std::vector<std::vector<double> > l(2,std::vector<double> (5));
        l[0][0]=1.0;
        l[1][0]=1.0;
        l[0][1]=2.0*qPos[0]-1.0;
        l[1][1]=2.0*qPos[1]-1.0;
        l[0][2]=6.0*qPos[0]*qPos[0]-6.0*qPos[0]+1.0;
        l[1][2]=6.0*qPos[1]*qPos[1]-6.0*qPos[1]+1.0;
        l[0][3]=20.0*qPos[0]*qPos[0]*qPos[0] - 30.0*qPos[0]*qPos[0] + 12.0*qPos[0] - 1.0;
        l[1][3]=20.0*qPos[1]*qPos[1]*qPos[1] - 30.0*qPos[1]*qPos[1] + 12.0*qPos[1] - 1.0;
        l[0][4]=1.0-20.0*qPos[0]+90.0*pow(qPos[0],2)-140.0*pow(qPos[0],3)+70.0*pow(qPos[0],4);
        l[1][4]=1.0-20.0*qPos[1]+90.0*pow(qPos[1],2)-140.0*pow(qPos[1],3)+70.0*pow(qPos[1],4);

        for (int i=0;i<4;i++)
          for (int j=0;j<5;j++)
            out[20+i*5+j] +=y[0]*l[0][i]*l[1][j]*it->weight();

        for (int i=0;i<5;i++)
          for (int j=0;j<4;j++)
            out[40+i*4+j] +=y[1]*l[0][i]*l[1][j]*it->weight();
      }
    }

  private:
    typename LB::Traits::RangeFieldType sign0, sign1, sign2, sign3;
    typename LB::Traits::DomainType n0, n1, n2, n3;
  };
}

#endif // DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS3_CUBE2D_LOCALINTERPOLATION_HH
