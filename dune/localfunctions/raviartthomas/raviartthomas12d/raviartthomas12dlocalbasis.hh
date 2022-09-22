// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS12DLOCALBASIS_HH
#define DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS12DLOCALBASIS_HH

#include <numeric>
#include <vector>

#include <dune/common/fmatrix.hh>

#include "../../common/localbasis.hh"

namespace Dune
{
  /**
   * @ingroup LocalBasisImplementation
   * \brief First order Raviart-Thomas shape functions on the reference triangle.
   *
   * \tparam D Type to represent the field in the domain.
   * \tparam R Type to represent the field in the range.
   *
   * \nosubgrouping
   */
  template<class D, class R>
  class RT12DLocalBasis
  {

  public:
    typedef LocalBasisTraits<D,2,Dune::FieldVector<D,2>,R,2,Dune::FieldVector<R,2>,
        Dune::FieldMatrix<R,2,2> > Traits;

    /**
     * \brief Make set number s, where 0 <= s < 8
     *
     * \param s Edge orientation indicator
     */
    RT12DLocalBasis (std::bitset<3> s = 0)
    {
      for (size_t i=0; i<3; i++)
        sign_[i] = (s[i]) ? -1.0 : 1.0;
    }

    //! \brief number of shape functions
    unsigned int size () const
    {
      return 8;
    }

    /**
     * \brief Evaluate all shape functions
     *
     * \param in Position
     * \param out return value
     */
    inline void evaluateFunction (const typename Traits::DomainType& in,
                                  std::vector<typename Traits::RangeType>& out) const
    {
      out.resize(8);
      out[0][0] = sign_[0]*(in[0] - 4.0*in[0]*in[1]);
      out[0][1] = sign_[0]*(-1.0 + 5.0*in[1] - 4.0*in[1]*in[1]);
      out[1][0] = sign_[1]*(-1.0 + 5.0*in[0] - 4.0*in[0]*in[0]);
      out[1][1] = sign_[1]*(in[1] - 4.0*in[0]*in[1]);
      out[2][0] = sign_[2]*(-3.0*in[0] + 4.0*in[0]*in[0] + 4.0*in[1]*in[0]);
      out[2][1] = sign_[2]*(-3.0*in[1] + 4.0*in[0]*in[1] + 4.0*in[1]*in[1]);
      out[3][0] = -5.0*in[0] + 8.0*in[0]*in[0] + 4.0*in[1]*in[0];
      out[3][1] = 3.0 - 6.0*in[0] - 7.0*in[1] + 8.0*in[0]*in[1] + 4.0*in[1]*in[1];
      out[4][0] = -3.0 + 7.0*in[0] + 6.0*in[1] - 4.0*in[0]*in[0] - 8.0*in[1]*in[0];
      out[4][1] = 5.0*in[1] - 4.0*in[0]*in[1] - 8.0*in[1]*in[1];
      out[5][0] = in[0] - 4.0*in[0]*in[0] + 4.0*in[1]*in[0];
      out[5][1] = -1.0*in[1] - 4.0*in[0]*in[1] + 4.0*in[1]*in[1];
      out[6][0] = 16.0*in[0] - 16.0*in[0]*in[0] - 8.0*in[1]*in[0];
      out[6][1] = 8.0*in[1] - 16.0*in[0]*in[1] - 8.0*in[1]*in[1];
      out[7][0] = 8.0*in[0] - 8.0*in[0]*in[0] - 16.0*in[1]*in[0];
      out[7][1] = 16.0*in[1] - 8.0*in[0]*in[1] - 16.0*in[1]*in[1];
    }

    /**
     * \brief Evaluate Jacobian of all shape functions
     *
     * \param in Position
     * \param out return value
     */
    inline void evaluateJacobian (const typename Traits::DomainType& in,
                                  std::vector<typename Traits::JacobianType>& out) const
    {
      out.resize(8);

      out[0][0][0] = sign_[0]*(1.0 - 4.0*in[1]);
      out[0][0][1] = sign_[0]*(-4.0*in[0]);
      out[0][1][0] = 0.0;
      out[0][1][1] = sign_[0]*(5.0 - 8.0*in[1]);

      out[1][0][0] = sign_[1]*(5.0 - 8.0*in[0]);
      out[1][0][1] = 0.0;
      out[1][1][0] = sign_[1]*(-4.0*in[1]);
      out[1][1][1] = sign_[1]*(1.0 - 4.0*in[0]);

      out[2][0][0] = sign_[2]*(-3.0 + 8.0*in[0] + 4.0*in[1]);
      out[2][0][1] = sign_[2]*(4.0*in[0]);
      out[2][1][0] = sign_[2]*(4.0*in[1]);
      out[2][1][1] = sign_[2]*(-3.0 + 4.0*in[0] + 8.0*in[1]);

      out[3][0][0] = -5.0 + 16.0*in[0] + 4.0*in[1];
      out[3][0][1] = 4.0*in[0];
      out[3][1][0] = -6.0 + 8.0*in[1];
      out[3][1][1] = -7.0 + 8.0*in[0] + 8.0*in[1];

      out[4][0][0] = 7.0 - 8.0*in[0] - 8.0*in[1];
      out[4][0][1] = 6.0 - 8.0*in[0];
      out[4][1][0] = -4.0*in[1];
      out[4][1][1] = 5.0 - 4.0*in[0] - 16.0*in[1];

      out[5][0][0] = 1.0 - 8.0*in[0] + 4*in[1];
      out[5][0][1] = 4.0*in[0];
      out[5][1][0] = -4.0*in[1];
      out[5][1][1] = -1.0 - 4.0*in[0] + 8.0*in[1];

      out[6][0][0] = 16.0 - 32.0*in[0] - 8.0*in[1];
      out[6][0][1] = -8.0*in[0];
      out[6][1][0] = -16.0*in[1];
      out[6][1][1] = 8.0 - 16.0*in[0] - 16.0*in[1];

      out[7][0][0] = 8.0 - 16.0*in[0] - 16.0*in[1];
      out[7][0][1] = -16.0*in[0];
      out[7][1][0] = -8.0*in[1];
      out[7][1][1] = 16.0 - 8.0*in[0] - 32.0*in[1];
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
        auto const direction = std::distance(order.begin(), std::find(order.begin(), order.end(), 1));
        out.resize(size());

        switch (direction) {
          case 0:
            out[0][0] = sign_[0]*(1.0 - 4.0*in[1]);
            out[0][1] = 0.0;
            out[1][0] = sign_[1]*(5.0 - 8.0*in[0]);
            out[1][1] = sign_[1]*(-4.0*in[1]);
            out[2][0] = sign_[2]*(-3.0 + 8.0*in[0] + 4.0*in[1]);
            out[2][1] = sign_[2]*(4.0*in[1]);
            out[3][0] = -5.0 + 16.0*in[0] + 4.0*in[1];
            out[3][1] = -6.0 + 8.0*in[1];
            out[4][0] = 7.0 - 8.0*in[0] - 8.0*in[1];
            out[4][1] = -4.0*in[1];
            out[5][0] = 1.0 - 8.0*in[0] + 4*in[1];
            out[5][1] = -4.0*in[1];
            out[6][0] = 16.0 - 32.0*in[0] - 8.0*in[1];
            out[6][1] = -16.0*in[1];
            out[7][0] = 8.0 - 16.0*in[0] - 16.0*in[1];
            out[7][1] = -8.0*in[1];
            break;
          case 1:
            out[2][1] = sign_[2]*(-3.0 + 4.0*in[0] + 8.0*in[1]);
            out[2][0] = sign_[2]*(4.0*in[0]);
            out[1][1] = sign_[1]*(1.0 - 4.0*in[0]);
            out[1][0] = 0.0;
            out[0][0] = sign_[0]*(-4.0*in[0]);
            out[0][1] = sign_[0]*(5.0 - 8.0*in[1]);
            out[3][0] = 4.0*in[0];
            out[3][1] = -7.0 + 8.0*in[0] + 8.0*in[1];
            out[4][0] = 6.0 - 8.0*in[0];
            out[4][1] = 5.0 - 4.0*in[0] - 16.0*in[1];
            out[5][0] = 4.0*in[0];
            out[5][1] = -1.0 - 4.0*in[0] + 8.0*in[1];
            out[6][0] = -8.0*in[0];
            out[6][1] = 8.0 - 16.0*in[0] - 16.0*in[1];
            out[7][0] = -16.0*in[0];
            out[7][1] = 16.0 - 8.0*in[0] - 32.0*in[1];
            break;
          default:
            DUNE_THROW(RangeError, "Component out of range.");
        }
      } else {
        DUNE_THROW(NotImplemented, "Desired derivative order is not implemented");
      }
    }

    //! \brief Polynomial order of the shape functions
    unsigned int order () const
    {
      return 2;
    }

  private:
    std::array<R,3> sign_;
  };
}
#endif // DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS12DLOCALBASIS_HH
