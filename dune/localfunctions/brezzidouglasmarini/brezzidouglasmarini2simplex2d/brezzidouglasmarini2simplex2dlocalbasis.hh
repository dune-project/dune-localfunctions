// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_BREZZIDOUGLASMARINI2_SIMPLEX2D_LOCALBASIS_HH
#define DUNE_LOCALFUNCTIONS_BREZZIDOUGLASMARINI2_SIMPLEX2D_LOCALBASIS_HH

#include <array>
#include <bitset>
#include <numeric>
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
  class BDM2Simplex2DLocalBasis
  {

  public:
    typedef LocalBasisTraits<D,2,Dune::FieldVector<D,2>,
        R,2,Dune::FieldVector<R,2>,
        Dune::FieldMatrix<R,2,2> > Traits;

    //! \brief Standard constructor
    BDM2Simplex2DLocalBasis()
    {
      for (size_t i=0; i<3; i++)
        sign_[i] = 1.0;
    }

    /**
     * \brief Make set number s, where 0 <= s < 8
     *
     * \param s Edge orientation indicator
     */
    BDM2Simplex2DLocalBasis(std::bitset<3> s)
    {
      for (size_t i=0; i<3; i++)
        sign_[i] = s[i] ? -1.0 : 1.0;
    }

    //! \brief number of shape functions
    unsigned int size() const
    {
      return 12;
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

      out[0][0] = sign_[0]*(-2*in[0]*in[1] + in[0]*in[0]);
      out[0][1] = sign_[0]*(-1 + 6*in[1] -2*in[0]*in[1] - 5*in[1]*in[1]);

      out[1][0] = 1.5*in[0] + 3*in[0]*in[1] - 4.5*in[0]*in[0];
      out[1][1] = -3 + 6*in[0] + 10.5*in[1] - 15*in[0]*in[1] - 7.5*in[1]*in[1];

      out[2][0] = sign_[0]*(-7.5*in[0] + 5*in[0]*in[1] + 12.5*in[0]*in[0]);
      out[2][1] = sign_[0]*(-5 + 30*in[0] + 7.5*in[1] - 25*in[0]*in[1] - 30*in[0]*in[0] - 2.5*in[1]*in[1]);



      out[3][0] = sign_[1]*(-1 + 6*in[0] - 2*in[0]*in[1] - 5*in[0]*in[0]);
      out[3][1] = sign_[1]*(-2*in[0]*in[1] + in[1]*in[1]);

      out[4][0] = 3 - 10.5*in[0] - 6*in[1] + 15*in[0]*in[1] + 7.5*in[0]*in[0];
      out[4][1] = -1.5*in[1] - 3*in[0]*in[1] + 4.5*in[1]*in[1];

      out[5][0] = sign_[1]*(-5 + 7.5*in[0] + 30*in[1] - 25*in[0]*in[1] - 2.5*in[0]*in[0] - 30*in[1]*in[1]);
      out[5][1] = sign_[1]*(-7.5*in[1] + 5*in[0]*in[1] + 12.5*in[1]*in[1]);



      out[6][0] = sign_[2]*(-3*in[0] + 4*in[0]*in[1] + 4*in[0]*in[0]);
      out[6][1] = sign_[2]*(-3*in[1] + 4*in[0]*in[1] + 4*in[1]*in[1]);

      out[7][0] = -3*in[0] + 6*in[0]*in[0];
      out[7][1] = 3*in[1] - 6*in[1]*in[1];

      out[8][0] = sign_[2]*(-10*in[0]*in[1] + 5*in[0]*in[0]);
      out[8][1] = sign_[2]*(-10*in[0]*in[1] + 5*in[1]*in[1]);



      out[9][0] = 18*in[0] - 12*in[0]*in[1] - 18*in[0]*in[0];
      out[9][1] = 6*in[1] - 12*in[0]*in[1] - 6*in[1]*in[1];

      out[10][0] = 6*in[0] - 12*in[0]*in[1] - 6*in[0]*in[0];
      out[10][1] = 18*in[1] - 12*in[0]*in[1] - 18*in[1]*in[1];

      out[11][0] = 90*in[0] - 180*in[0]*in[1] - 90*in[0]*in[0];
      out[11][1] = -90*in[1] + 180*in[0]*in[1] + 90*in[1]*in[1];
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

      out[0][0][0] = sign_[0]*(-2*in[1] + 2*in[0]);
      out[0][0][1] = sign_[0]*(-2*in[0]);

      out[0][1][0] = sign_[0]*(-2*in[1]);
      out[0][1][1] = sign_[0]*(6 -2*in[0] - 10*in[1]);


      out[1][0][0] = 1.5 + 3*in[1] - 9*in[0];
      out[1][0][1] = 3*in[0];

      out[1][1][0] = 6 - 15*in[1];
      out[1][1][1] = 10.5 - 15*in[0] - 15*in[1];


      out[2][0][0] = sign_[0]*(-7.5 + 5*in[1] + 25*in[0]);
      out[2][0][1] = sign_[0]*(5*in[0]);

      out[2][1][0] = sign_[0]*(30 - 25*in[1] - 60*in[0]);
      out[2][1][1] = sign_[0]*(7.5 - 25*in[0] - 5*in[1]);



      out[3][0][0] = sign_[1]*(6 - 2*in[1] - 10*in[0]);
      out[3][0][1] = sign_[1]*(-2*in[0]);

      out[3][1][0] = sign_[1]*(-2*in[1]);
      out[3][1][1] = sign_[1]*(-2*in[0] + 2*in[1]);


      out[4][0][0] = -10.5 + 15*in[1] + 15*in[0];
      out[4][0][1] = -6 + 15*in[0];

      out[4][1][0] = -3*in[1];
      out[4][1][1] = -1.5 - 3*in[0] + 9*in[1];


      out[5][0][0] = sign_[1]*(7.5 - 25*in[1] - 5*in[0]);
      out[5][0][1] = sign_[1]*(30 - 25*in[0] - 60*in[1]);

      out[5][1][0] = sign_[1]*(5*in[1]);
      out[5][1][1] = sign_[1]*(-7.5 + 5*in[0] + 25*in[1]);



      out[6][0][0] = sign_[2]*(-3 + 4*in[1] + 8*in[0]);
      out[6][0][1] = sign_[2]*(4*in[0]);

      out[6][1][0] = sign_[2]*(4*in[1]);
      out[6][1][1] = sign_[2]*(-3 + 4*in[0] + 8*in[1]);


      out[7][0][0] = -3 + 12*in[0];
      out[7][0][1] = 0;

      out[7][1][0] = 0;
      out[7][1][1] = 3 - 12*in[1];


      out[8][0][0] = sign_[2]*(-10*in[1] + 10*in[0]);
      out[8][0][1] = sign_[2]*(-10*in[0]);

      out[8][1][0] = sign_[2]*(-10*in[1]);
      out[8][1][1] = sign_[2]*(-10*in[0] + 10*in[1]);


      out[9][0][0] = 18 - 12*in[1] - 36*in[0];
      out[9][0][1] = -12*in[0];

      out[9][1][0] = -12*in[1];
      out[9][1][1] = 6 - 12*in[0] - 12*in[1];

      out[10][0][0] = 6 - 12*in[1] - 12*in[0];
      out[10][0][1] = -12*in[0];

      out[10][1][0] = -12*in[1];
      out[10][1][1] = 18 - 12*in[0] - 36*in[1];

      out[11][0][0] = 90 - 180*in[1] - 180*in[0];
      out[11][0][1] = -180*in[0];

      out[11][1][0] = 180*in[1];
      out[11][1][1] = -90 + 180*in[0] + 180*in[1];
    }

    //! \brief Evaluate partial derivatives of all shape functions
    void partial (const std::array<unsigned int, 2>& order,
                  const typename Traits::DomainType& in,         // position
                  std::vector<typename Traits::RangeType>& out) const      // return value
    {
      auto totalOrder = std::accumulate(order.begin(), order.end(), 0);
      if (totalOrder == 0) {
        evaluateFunction(in, out);
      } else if (totalOrder == 1) {
        out.resize(size());
        auto const direction = std::distance(order.begin(), std::find(order.begin(), order.end(), 1));

        switch (direction) {
        case 0:
          out[0][0] = sign_[0]*(-2*in[1] + 2*in[0]);
          out[0][1] = sign_[0]*(-2*in[1]);

          out[1][0] = 1.5 + 3*in[1] - 9*in[0];
          out[1][1] = 6 - 15*in[1];

          out[2][0] = sign_[0]*(-7.5 + 5*in[1] + 25*in[0]);
          out[2][1] = sign_[0]*(30 - 25*in[1] - 60*in[0]);

          out[3][0] = sign_[1]*(6 - 2*in[1] - 10*in[0]);
          out[3][1] = sign_[1]*(-2*in[1]);

          out[4][0] = -10.5 + 15*in[1] + 15*in[0];
          out[4][1] = -3*in[1];

          out[5][0] = sign_[1]*(7.5 - 25*in[1] - 5*in[0]);
          out[5][1] = sign_[1]*(5*in[1]);

          out[6][0] = sign_[2]*(-3 + 4*in[1] + 8*in[0]);
          out[6][1] = sign_[2]*(4*in[1]);

          out[7][0] = -3 + 12*in[0];
          out[7][1] = 0;

          out[8][0] = sign_[2]*(-10*in[1] + 10*in[0]);
          out[8][1] = sign_[2]*(-10*in[1]);

          out[9][0] = 18 - 12*in[1] - 36*in[0];
          out[9][1] = -12*in[1];

          out[10][0] = 6 - 12*in[1] - 12*in[0];
          out[10][1] = -12*in[1];

          out[11][0] = 90 - 180*in[1] - 180*in[0];
          out[11][1] = 180*in[1];
          break;
        case 1:
          out[0][0] = sign_[0]*(-2*in[0]);
          out[0][1] = sign_[0]*(6 -2*in[0] - 10*in[1]);

          out[1][0] = 3*in[0];
          out[1][1] = 10.5 - 15*in[0] - 15*in[1];

          out[2][0] = sign_[0]*(5*in[0]);
          out[2][1] = sign_[0]*(7.5 - 25*in[0] - 5*in[1]);

          out[3][0] = sign_[1]*(-2*in[0]);
          out[3][1] = sign_[1]*(-2*in[0] + 2*in[1]);

          out[4][0] = -6 + 15*in[0];
          out[4][1] = -1.5 - 3*in[0] + 9*in[1];

          out[5][0] = sign_[1]*(30 - 25*in[0] - 60*in[1]);
          out[5][1] = sign_[1]*(-7.5 + 5*in[0] + 25*in[1]);

          out[6][0] = sign_[2]*(4*in[0]);
          out[6][1] = sign_[2]*(-3 + 4*in[0] + 8*in[1]);

          out[7][0] = 0;
          out[7][1] = 3 - 12*in[1];

          out[8][0] = sign_[2]*(-10*in[0]);
          out[8][1] = sign_[2]*(-10*in[0] + 10*in[1]);

          out[9][0] = -12*in[0];
          out[9][1] = 6 - 12*in[0] - 12*in[1];

          out[10][0] = -12*in[0];
          out[10][1] = 18 - 12*in[0] - 36*in[1];

          out[11][0] = -180*in[0];
          out[11][1] = -90 + 180*in[0] + 180*in[1];
          break;
        default:
          DUNE_THROW(RangeError, "Component out of range.");
        }
      } else {
        DUNE_THROW(NotImplemented, "Desired derivative order is not implemented");
      }
    }

    //! \brief Polynomial order of the shape functions
    unsigned int order() const
    {
      return 2; // TODO: check whether this is not order 3
    }

  private:
    std::array<R,3> sign_;
  };
} // end namespace Dune
#endif // DUNE_LOCALFUNCTIONS_BREZZIDOUGLASMARINI2_SIMPLEX2D_LOCALBASIS_HH
