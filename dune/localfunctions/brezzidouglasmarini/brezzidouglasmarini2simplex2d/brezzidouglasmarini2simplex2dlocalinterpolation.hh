// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LOCALFUNCTIONS_BREZZIDOUGLASMARINI2_SIMPLEX2D_LOCALINTERPOLATION_HH
#define DUNE_LOCALFUNCTIONS_BREZZIDOUGLASMARINI2_SIMPLEX2D_LOCALINTERPOLATION_HH

#include <vector>

#include <dune/geometry/quadraturerules.hh>

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
     * \param f function which should be interpolated
     * \param out return value, vector of coefficients
     */
    template<typename F, typename C>
    void interpolate(const F& f, std::vector<C>& out) const
    {
      // f gives v*outer normal at a point on the edge!
      typedef typename LB::Traits::RangeFieldType Scalar;
      typename F::Traits::RangeType y;

      out.resize(12);
      fill(out.begin(), out.end(), 0.0);

      const int qOrder = 4;
      const Dune::QuadratureRule<Scalar,1>& rule = Dune::QuadratureRules<Scalar,1>::rule(Dune::GeometryType(Dune::GeometryType::simplex,1), qOrder);

      for (typename Dune::QuadratureRule<Scalar,1>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
      {
        // TODO: write interpolation
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
