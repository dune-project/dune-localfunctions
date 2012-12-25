// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS1_CUBE3D_LOCALINTERPOLATION_HH
#define DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS1_CUBE3D_LOCALINTERPOLATION_HH

#include <vector>

#include <dune/geometry/quadraturerules.hh>

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
    //! \brief Standard constructor
    RT1Cube3DLocalInterpolation ()
    {
      sign0 = sign1 = sign2 = sign3 = sign4 = sign5 = 1.0;
    }

    /**
     * \brief Make set number s, where 0 <= s < 64
     *
     * \param s Edge orientation indicator
     */
    RT1Cube3DLocalInterpolation (unsigned int s)
    {
      sign0 = sign1 = sign2 = sign3 = sign4 = sign5 = 1.0;
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
      if (s & 8)
      {
        sign3 = -1.0;
      }
      if (s & 16)
      {
        sign4 = -1.0;
      }
      if (s & 32)
      {
        sign5 = -1.0;
      }

      n0[0] = -1.0;
      n0[1] =  0.0;
      n0[2] =  0.0;
      n1[0] =  1.0;
      n1[1] =  0.0;
      n1[2] =  0.0;
      n2[0] =  0.0;
      n2[1] = -1.0;
      n2[2] =  0.0;
      n3[0] =  0.0;
      n3[1] =  1.0;
      n3[2] =  0.0;
      n4[0] =  0.0;
      n4[1] =  0.0;
      n4[2] = -1.0;
      n5[0] =  0.0;
      n5[1] =  0.0;
      n5[2] =  1.0;
    }

    /**
     * \brief Interpolate a given function with shape functions
     *
     * \tparam F Function type for function which should be interpolated
     * \tparam C Coefficient type
     * \param f function which should be interpolated
     * \param out return value, vector of coefficients
     */
    template<class F, class C>
    void interpolate (const F& f, std::vector<C>& out) const
    {
      // f gives v*outer normal at a point on the edge!
      typedef typename LB::Traits::RangeFieldType Scalar;
      typedef typename LB::Traits::DomainFieldType Vector;
      typename F::Traits::RangeType y;

      out.resize(36);
      fill(out.begin(), out.end(), 0.0);

      const int qOrder = 3;
      const QuadratureRule<Scalar,2>& rule1 = QuadratureRules<Scalar,2>::rule(GeometryType(GeometryType::cube,2), qOrder);

      for (typename QuadratureRule<Scalar,2>::const_iterator it = rule1.begin();
           it != rule1.end(); ++it)
      {
        Dune::FieldVector<Scalar,2> qPos = it->position();
        typename LB::Traits::DomainType localPos;

        localPos[0] = 0.0;
        localPos[1] = qPos[0];
        localPos[2] = qPos[1];
        f.evaluate(localPos, y);
        out[0] += (y[0]*n0[0] + y[1]*n0[1] + y[2]*n0[2])*it->weight()*sign0;
        out[6] += (y[0]*n0[0] + y[1]*n0[1] + y[2]*n0[2])*(2.0*qPos[0] - 1.0)*it->weight();
        out[12] += (y[0]*n0[0] + y[1]*n0[1] + y[2]*n0[2])*(2.0*qPos[1] - 1.0)*it->weight();
        out[18] += (y[0]*n0[0] + y[1]*n0[1] + y[2]*n0[2])*(2.0*qPos[0] - 1.0)*(2.0*qPos[1] - 1.0)*it->weight();

        localPos[0] = 1.0;
        localPos[1] = qPos[0];
        localPos[2] = qPos[1];
        f.evaluate(localPos, y);
        out[1] += (y[0]*n1[0] + y[1]*n1[1] + y[2]*n1[2])*it->weight()*sign1;
        out[7] += (y[0]*n1[0] + y[1]*n1[1] + y[2]*n1[2])*(1.0 - 2.0*qPos[0])*it->weight();
        out[13] += (y[0]*n1[0] + y[1]*n1[1] + y[2]*n1[2])*(1.0 - 2.0*qPos[1])*it->weight();
        out[19] += (y[0]*n1[0] + y[1]*n1[1] + y[2]*n1[2])*(1.0 - 2.0*qPos[0])*(2.0*qPos[1] - 1.0)*it->weight();

        localPos[0] = qPos[0];
        localPos[1] = 0.0;
        localPos[2] = qPos[1];
        f.evaluate(localPos, y);
        out[2] += (y[0]*n2[0] + y[1]*n2[1] + y[2]*n2[2])*it->weight()*sign2;
        out[8] += (y[0]*n2[0] + y[1]*n2[1] + y[2]*n2[2])*(1.0 - 2.0*qPos[0])*it->weight();
        out[14] += (y[0]*n2[0] + y[1]*n2[1] + y[2]*n2[2])*(2.0*qPos[1] - 1.0)*it->weight();
        out[20] += (y[0]*n2[0] + y[1]*n2[1] + y[2]*n2[2])*(1.0 - 2.0*qPos[0])*(2.0*qPos[1] - 1.0)*it->weight();

        localPos[0] = qPos[0];
        localPos[1] = 1.0;
        localPos[2] = qPos[1];
        f.evaluate(localPos, y);
        out[3] += (y[0]*n3[0] + y[1]*n3[1] + y[2]*n3[2])*it->weight()*sign3;
        out[9] += (y[0]*n3[0] + y[1]*n3[1] + y[2]*n3[2])*(2.0*qPos[0] - 1.0)*it->weight();
        out[15] += (y[0]*n3[0] + y[1]*n3[1] + y[2]*n3[2])*(1.0 - 2.0*qPos[1])*it->weight();
        out[21] += (y[0]*n3[0] + y[1]*n3[1] + y[2]*n3[2])*(2.0*qPos[0] - 1.0)*(2.0*qPos[1] - 1.0)*it->weight();

        localPos[0] = qPos[0];
        localPos[1] = qPos[1];
        localPos[2] = 0.0;
        f.evaluate(localPos, y);
        out[4] += (y[0]*n4[0] + y[1]*n4[1] + y[2]*n4[2])*it->weight()*sign4;
        out[10] += (y[0]*n4[0] + y[1]*n4[1] + y[2]*n4[2])*(1.0 - 2.0*qPos[0])*it->weight();
        out[16] += (y[0]*n4[0] + y[1]*n4[1] + y[2]*n4[2])*(1.0 - 2.0*qPos[1])*it->weight();
        out[22] += (y[0]*n4[0] + y[1]*n4[1] + y[2]*n4[2])*(1.0 - 2.0*qPos[0])*(2.0*qPos[1] - 1.0)*it->weight();

        localPos[0] = qPos[0];
        localPos[1] = qPos[1];
        localPos[2] = 1.0;
        f.evaluate(localPos, y);
        out[5] += (y[0]*n5[0] + y[1]*n5[1] + y[2]*n5[2])*it->weight()*sign5;
        out[11] += (y[0]*n5[0] + y[1]*n5[1] + y[2]*n5[2])*(2.0*qPos[0] - 1.0)*it->weight();
        out[17] += (y[0]*n5[0] + y[1]*n5[1] + y[2]*n5[2])*(2.0*qPos[1] - 1.0)*it->weight();
        out[23] += (y[0]*n5[0] + y[1]*n5[1] + y[2]*n5[2])*(2.0*qPos[0] - 1.0)*(2.0*qPos[1] - 1.0)*it->weight();
      }

      const QuadratureRule<Vector,3>& rule2 = QuadratureRules<Vector,3>::rule(GeometryType(GeometryType::cube,3), qOrder);
      for (typename QuadratureRule<Vector,3>::const_iterator it = rule2.begin();
           it != rule2.end(); ++it)
      {
        FieldVector<double,3> qPos = it->position();

        f.evaluate(qPos, y);
        out[24] += y[0]*it->weight();
        out[25] += y[1]*it->weight();
        out[26] += y[2]*it->weight();
        out[27] += y[0]*qPos[1]*it->weight();
        out[28] += y[0]*qPos[2]*it->weight();
        out[29] += y[1]*qPos[0]*it->weight();
        out[30] += y[1]*qPos[2]*it->weight();
        out[31] += y[2]*qPos[0]*it->weight();
        out[32] += y[2]*qPos[1]*it->weight();
        out[33] += y[0]*qPos[1]*qPos[2]*it->weight();
        out[34] += y[1]*qPos[0]*qPos[2]*it->weight();
        out[35] += y[2]*qPos[0]*qPos[1]*it->weight();
      }
    }

  private:
    typename LB::Traits::RangeFieldType sign0, sign1, sign2, sign3, sign4, sign5;
    typename LB::Traits::DomainType n0, n1, n2, n3, n4, n5;
  };
}
#endif // DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS1_CUBE3D_LOCALINTERPOLATION_HH
