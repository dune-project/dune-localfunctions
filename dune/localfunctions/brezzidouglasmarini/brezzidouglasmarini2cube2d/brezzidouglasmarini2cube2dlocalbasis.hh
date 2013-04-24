// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LOCALFUNCTIONS_BREZZIDOUGLASMARINI2_CUBE2D_LOCALBASIS_HH
#define DUNE_LOCALFUNCTIONS_BREZZIDOUGLASMARINI2_CUBE2D_LOCALBASIS_HH

#include <vector>

#include <dune/common/fmatrix.hh>

#include "../../common/localbasis.hh"

namespace Dune
{
  /**
   * \ingroup LocalBasisImplementation
   * \brief First order Brezzi-Douglas-Marini shape functions on quadrilaterals.
   *
   * \tparam D Type to represent the field in the domain.
   * \tparam R Type to represent the field in the range.
   *
   * \nosubgrouping
   */
  template<class D, class R>
  class BDM2Cube2DLocalBasis
  {

  public:
    typedef LocalBasisTraits<D,2,Dune::FieldVector<D,2>,
        R,2,Dune::FieldVector<R,2>,
        Dune::FieldMatrix<R,2,2> > Traits;

    //! \brief Standard constructor
    BDM2Cube2DLocalBasis()
    {
      sign0 = sign1 = sign2 = sign3 = 1.0;
    }

    /**
     * \brief Make set number s, where 0 <= s < 16
     *
     * \param s Edge orientation indicator
     */
    BDM2Cube2DLocalBasis(unsigned int s)
    {
      sign0 = sign1 = sign2 = sign3 = 1.0;
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
    }

    //! \brief number of shape functions
    unsigned int size() const
    {
      return 14;
    }

    /**
     * \brief Evaluate all shape functions
     *
     * \param in Position
     * \param out return value
     */
    inline void evaluateFunction(const typename Traits::DomainType& in,
                                 std::vector<typename Traits::RangeType>& out) const
    {
      out.resize(size());

      out[0][0] = sign0*(-2.25 + 5.25*in[0] + 7.5*in[1] - 7.5*in[0]*in[1] - 3.0*in[0]*in[0] - 7.5*in[1]*in[1] + 7.5*in[0]*in[1]*in[1]);
      out[0][1] = sign0*(-1.25*in[1] + 3.75*in[1]*in[1] - 2.5*in[1]*in[1]*in[1]);
      out[1][0] = 3.0 - 3.0*in[0]-6.0*in[1] + 6.0*in[0]*in[1];
      out[1][1] = 0.0;
      out[2][0] = sign0*(-3.75 + 3.75*in[0] + 22.5*in[1] - 22.5*in[0]*in[1] - 22.5*in[1]*in[1] + 22.5*in[0]*in[1]*in[1]);
      out[2][1] = sign0*(-3.75*in[1] + 11.25*in[1]*in[1] - 7.5*in[1]*in[1]*in[1]);
      out[3][0] = sign1*(-0.75*in[0] - 7.5*in[0]*in[1] + 3.0*in[0]*in[0] + 7.5*in[0]*in[1]*in[1]);
      out[3][1] = sign1*(-1.25*in[1] + 3.75*in[1]*in[1] - 2.5*in[1]*in[1]*in[1]);
      out[4][0] = 3.0*in[0] - 6.0*in[0]*in[1];
      out[4][1] = 0.0;
      out[5][0] = sign1*(+3.75*in[0] - 22.5*in[0]*in[1] + 22.5*in[0]*in[1]*in[1]);
      out[5][1] = sign1*(-3.75*in[1] + 11.25*in[1]*in[1] - 7.5*in[1]*in[1]*in[1]);
      out[6][0] = sign2*(-1.25*in[0] + 3.75*in[0]*in[0] - 2.5*in[0]*in[0]*in[0]);
      out[6][1] = sign2*(-2.25 + 7.5*in[0] + 5.25*in[1] - 7.5*in[0]*in[1] - 7.5*in[0]*in[0] - 3.0*in[1]*in[1] + 7.5*in[0]*in[0]*in[1]);
      out[7][0] = 0.0;
      out[7][1] = -3.0 + 6.0*in[0] + 3.0*in[1] - 6.0*in[0]*in[1];
      out[8][0] = sign2*(-3.75*in[0] + 11.25*in[0]*in[0] - 7.5*in[0]*in[0]*in[0]);
      out[8][1] = sign2*(-3.75 + 22.5*in[0] + 3.75*in[1] - 22.5*in[0]*in[1] - 22.5*in[0]*in[0] + 22.5*in[0]*in[0]*in[1]);
      out[9][0] = sign3*(-1.25*in[0] + 3.75*in[0]*in[0] - 2.5*in[0]*in[0]*in[0]);
      out[9][1] = sign3*(-0.75*in[1] - 7.5*in[0]*in[1] + 3.0*in[1]*in[1] + 7.5*in[0]*in[0]*in[1]);
      out[10][0] = 0.0;
      out[10][1] = -3.0*in[1] + 6.0*in[0]*in[1];
      out[11][0] = sign3*(-3.75*in[0] + 11.25*in[0]*in[0] - 7.5*in[0]*in[0]*in[0]);
      out[11][1] = sign3*(3.75*in[1] - 22.5*in[0]*in[1] + 22.5*in[0]*in[0]*in[1]);
      out[12][0] = 6.0*in[0] - 6.0*in[0]*in[0];
      out[12][1] = 0.0;
      out[13][0] = 0.0;
      out[13][1] = 6.0*in[1] - 6.0*in[1]*in[1];
    }

    /**
     * \brief Evaluate Jacobian of all shape functions
     *
     * \param in Position
     * \param out return value
     */
    inline void evaluateJacobian(const typename Traits::DomainType& in,
                                 std::vector<typename Traits::JacobianType>& out) const
    {
      out.resize(size());

      out[0][0][0] = sign0*(5.25 - 7.5*in[1] - 6.0*in[0] + 7.5*in[1]*in[1]);
      out[0][0][1] = sign0*(7.5 - 7.5*in[0] - 15.0*in[1] + 15.0*in[0]*in[1]);
      out[0][1][0] = 0.0;
      out[0][1][1] = sign0*(-1.25 + 7.5*in[1] - 7.5*in[1]*in[1]);

      out[1][0][0] = -3.0 + 6.0*in[1];
      out[1][0][1] = -6.0 + 6.0*in[0];
      out[1][1][0] = 0.0;
      out[1][1][1] = 0.0;

      out[2][0][0] = sign0*(3.75 - 22.5*in[1] + 22.5*in[1]*in[1]);
      out[2][0][1] = sign0*(22.5 - 22.5*in[0] - 45.0*in[1] + 45.0*in[0]*in[1]);
      out[2][1][0] = 0.0;
      out[2][1][1] = sign0*(-3.75 + 22.5*in[1] - 22.5*in[1]*in[1]);

      out[3][0][0] = sign1*(-0.75 - 7.5*in[1] + 6.0*in[0] + 7.5*in[1]*in[1]);
      out[3][0][1] = sign1*(-7.5*in[0] + 15.0*in[0]*in[1]);
      out[3][1][0] = 0.0;
      out[3][1][1] = sign1*(-1.25 + 7.5*in[1] - 7.5*in[1]*in[1]);

      out[4][0][0] = 3.0 - 6.0*in[1];
      out[4][0][1] = -6.0*in[0];
      out[4][1][0] = 0.0;
      out[4][1][1] = 0.0;

      out[5][0][0] = sign1*(3.75 - 22.5*in[1] + 22.5*in[1]*in[1]);
      out[5][0][1] = sign1*(-22.5*in[0] + 45.0*in[0]*in[1]);
      out[5][1][0] = 0.0;
      out[5][1][1] = sign1*(-3.75 + 22.5*in[1] - 22.5*in[1]*in[1]);

      out[6][0][0] = sign2*(-1.25 + 7.5*in[0] - 7.5*in[0]*in[0]);
      out[6][0][1] = 0.0;
      out[6][1][0] = sign2*(7.5 - 7.5*in[1] - 15.0*in[0] + 15.0*in[0]*in[1]);
      out[6][1][1] = sign2*(5.25 - 7.5*in[0]- 6.0*in[1] + 7.5*in[0]*in[0]);

      out[7][0][0] = 0.0;
      out[7][0][1] = 0.0;
      out[7][1][0] = 6.0 - 6.0*in[1];
      out[7][1][1] = 3.0 - 6.0*in[0];

      out[8][0][0] = sign2*(-3.75 + 22.5*in[0] - 22.5*in[0]*in[0]);
      out[8][0][1] = 0.0;
      out[8][1][0] = sign2*(22.5 - 22.5*in[1] - 45.0*in[0] + 45.0*in[0]*in[1]);
      out[8][1][1] = sign2*(3.75 - 22.5*in[0] + 22.5*in[0]*in[0]);

      out[9][0][0] = sign3*(-1.25 + 7.5*in[0] - 7.5*in[0]*in[0]);
      out[9][0][1] = 0.0;
      out[9][1][0] = sign3*(-7.5*in[1] + 15.0*in[0]*in[1]);
      out[9][1][1] = sign3*(-0.75 - 7.5*in[0] + 6.0*in[1] + 7.5*in[0]*in[0]);

      out[10][0][0] = 0.0;
      out[10][0][1] = 0.0;
      out[10][1][0] = 6.0*in[1];
      out[10][1][1] = -3.0 + 6.0*in[0];

      out[11][0][0] = sign3*(-3.75 + 22.5*in[0] - 22.5*in[0]*in[0]);
      out[11][0][1] = 0.0;
      out[11][1][0] = sign3*(-22.5*in[1] + 45*in[0]*in[1]);
      out[11][1][1] = sign3*(3.75 - 22.5*in[0] + 22.5*in[0]*in[0]);

      out[12][0][0] = 6.0 - 12.0*in[0];
      out[12][0][1] = 0.0;
      out[12][1][0] = 0.0;
      out[12][1][1] = 0.0;

      out[13][0][0] = 0.0;
      out[13][0][1] = 0.0;
      out[13][1][0] = 0.0;
      out[13][1][1] = 6.0 - 12.0*in[1];
    }

    //! \brief Polynomial order of the shape functions
    unsigned int order() const
    {
      return 3;
    }

  private:
    R sign0, sign1, sign2, sign3;
  };
} // end namespace Dune
#endif // DUNE_LOCALFUNCTIONS_BREZZIDOUGLASMARINI2_CUBE2D_LOCALBASIS_HH
