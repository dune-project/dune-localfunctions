// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS2_CUBE2D_LOCALINTERPOLATION_HH
#define DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS2_CUBE2D_LOCALINTERPOLATION_HH

#include <vector>

#include <dune/geometry/quadraturerules.hh>

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
    //! \brief Standard constructor
    RT2Cube2DLocalInterpolation ()
    {
      sign0 = sign1 = sign2 = sign3 = 1.0;
    }

    /**
     * \brief Make set number s, where 0 <= s < 8
     *
     * \param s Edge orientation indicator
     */
    RT2Cube2DLocalInterpolation (unsigned int s)
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
     * \param f function which should be interpolated
     * \param out return value, vector of coefficients
     */
    template<typename F, typename C>
    void interpolate (const F& f, std::vector<C>& out) const
    {
      // f gives v*outer normal at a point on the edge!
      typedef typename LB::Traits::RangeFieldType Scalar;
      typedef typename LB::Traits::DomainFieldType Vector;
      typename F::Traits::RangeType y;

      out.resize(24);
      fill(out.begin(), out.end(), 0.0);

      const int qOrder = 6;
      const QuadratureRule<Scalar,1>& rule = QuadratureRules<Scalar,1>::rule(GeometryType(GeometryType::cube,1), qOrder);

      for (typename QuadratureRule<Scalar,1>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
      {
        Scalar qPos = it->position();
        typename LB::Traits::DomainType localPos;

        localPos[0] = 0.0;
        localPos[1] = qPos;
        f.evaluate(localPos, y);
        out[0] += (y[0]*n0[0] + y[1]*n0[1])*it->weight()*sign0;
        out[1] += (y[0]*n0[0] + y[1]*n0[1])*(2.0*qPos - 1.0)*it->weight();
        out[2] += (y[0]*n0[0] + y[1]*n0[1])*(6.0*qPos*qPos - 6.0*qPos + 1.0)*it->weight()*sign0;

        localPos[0] = 1.0;
        localPos[1] = qPos;
        f.evaluate(localPos, y);
        out[3] += (y[0]*n1[0] + y[1]*n1[1])*it->weight()*sign1;
        out[4] += (y[0]*n1[0] + y[1]*n1[1])*(1.0 - 2.0*qPos)*it->weight();
        out[5] += (y[0]*n1[0] + y[1]*n1[1])*(6.0*qPos*qPos - 6.0*qPos + 1.0)*it->weight()*sign1;

        localPos[0] = qPos;
        localPos[1] = 0.0;
        f.evaluate(localPos, y);
        out[6] += (y[0]*n2[0] + y[1]*n2[1])*it->weight()*sign2;
        out[7] += (y[0]*n2[0] + y[1]*n2[1])*(1.0 - 2.0*qPos)*it->weight();
        out[8] += (y[0]*n2[0] + y[1]*n2[1])*(6.0*qPos*qPos - 6.0*qPos + 1.0)*it->weight()*sign2;

        localPos[0] = qPos;
        localPos[1] = 1.0;
        f.evaluate(localPos, y);
        out[9]  += (y[0]*n3[0] + y[1]*n3[1])*it->weight()*sign3;
        out[10] += (y[0]*n3[0] + y[1]*n3[1])*(2.0*qPos - 1.0)*it->weight();
        out[11] += (y[0]*n3[0] + y[1]*n3[1])*(6.0*qPos*qPos - 6.0*qPos + 1.0)*it->weight()*sign3;
      }

      const QuadratureRule<Vector,2>& rule2 = QuadratureRules<Vector,2>::rule(GeometryType(GeometryType::cube,2), qOrder);

      for (typename QuadratureRule<Vector,2>::const_iterator it = rule2.begin();
           it != rule2.end(); ++it)
      {
        FieldVector<double,2> qPos = it->position();

        f.evaluate(qPos, y);
        out[12] += y[0]*it->weight();
        out[13] += y[1]*it->weight();
        out[14] += y[0]*qPos[0]*it->weight();
        out[15] += y[1]*qPos[0]*it->weight();
        out[16] += y[0]*qPos[1]*it->weight();
        out[17] += y[1]*qPos[1]*it->weight();
        out[18] += y[0]*qPos[0]*qPos[1]*it->weight();
        out[19] += y[1]*qPos[0]*qPos[1]*it->weight();
        out[20] += y[0]*qPos[1]*qPos[1]*it->weight();
        out[21] += y[1]*qPos[0]*qPos[0]*it->weight();
        out[22] += y[0]*qPos[0]*qPos[1]*qPos[1]*it->weight();
        out[23] += y[1]*qPos[0]*qPos[0]*qPos[1]*it->weight();
      }
    }

  private:
    typename LB::Traits::RangeFieldType sign0, sign1, sign2, sign3;
    typename LB::Traits::DomainType n0, n1, n2, n3;
  };
}
#endif // DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS2_CUBE2D_LOCALINTERPOLATION_HH
