// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_BREZZIDOUGLASMARINI1_CUBE3D_LOCALINTERPOLATION_HH
#define DUNE_LOCALFUNCTIONS_BREZZIDOUGLASMARINI1_CUBE3D_LOCALINTERPOLATION_HH

#include <vector>

#include <dune/geometry/quadraturerules.hh>

namespace Dune
{

  /**
   * \brief First order Brezzi-Douglas-Marini shape functions on the reference
   *        hexahedron.
   *
   * \tparam LB corresponding LocalBasis giving traits
   *
   * \ingroup LocalInterpolationImplementation
   * \nosubgrouping
   */
  template<class LB>
  class BDM1Cube3DLocalInterpolation
  {

  public:
    //! \brief Standard constructor
    BDM1Cube3DLocalInterpolation()
    {
      sign0 = sign1 = sign2 = sign3 = sign4 = sign5 = 1.0;
    }

    /**
     * \brief Make set number s, where 0 <= s < 64
     *
     * \param s Edge orientation indicator
     */
    BDM1Cube3DLocalInterpolation(unsigned int s)
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
    template<typename F, typename C>
    void interpolate(const F& f, std::vector<C>& out) const
    {
      // f gives v*outer normal at a point on the edge!
      typedef typename LB::Traits::RangeFieldType Scalar;
      //typedef typename LB::Traits::DomainFieldType Vector;

      DUNE_THROW( NotImplemented, "Interpolation for BDM1Cube3D finite elements is not implemented." );

      out.resize(18);
      fill(out.begin(), out.end(), 0.0);

      const int qOrder = 4;
      const QuadratureRule<Scalar,1>& rule = QuadratureRules<Scalar,1>::rule(GeometryTypes::cube(1), qOrder);

      for (typename QuadratureRule<Scalar,1>::const_iterator it = rule.begin();
           it != rule.end(); ++it)
      {
        // TODO: write interpolation
      }
    }

  private:
    typename LB::Traits::RangeFieldType sign0, sign1, sign2, sign3, sign4, sign5;
    typename LB::Traits::DomainType n0, n1, n2, n3, n4, n5;
  };
} // end namespace Dune
#endif // DUNE_LOCALFUNCTIONS_BREZZIDOUGLASMARINI1_CUBE3D_LOCALINTERPOLATION_HH
