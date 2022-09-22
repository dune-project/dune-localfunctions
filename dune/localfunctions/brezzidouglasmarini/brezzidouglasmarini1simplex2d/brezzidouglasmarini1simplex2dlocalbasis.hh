// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_BREZZIDOUGLASMARINI1_SIMPLEX2D_LOCALBASIS_HH
#define DUNE_LOCALFUNCTIONS_BREZZIDOUGLASMARINI1_SIMPLEX2D_LOCALBASIS_HH

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
   * \brief First order Brezzi-Douglas-Marini shape functions on the reference triangle.
   *
   * \tparam D Type to represent the field in the domain.
   * \tparam R Type to represent the field in the range.
   *
   * \nosubgrouping
   */
  template<class D, class R>
  class BDM1Simplex2DLocalBasis
  {

  public:
    typedef LocalBasisTraits<D,2,Dune::FieldVector<D,2>,R,2,Dune::FieldVector<R,2>,
        Dune::FieldMatrix<R,2,2> > Traits;

    //! \brief Standard constructor
    BDM1Simplex2DLocalBasis ()
    {
      for (size_t i=0; i<3; i++)
        sign_[i] = 1.0;
    }

    /**
     * \brief Make set number s, where 0 <= s < 8
     *
     * \param s Edge orientation indicator
     */
    BDM1Simplex2DLocalBasis (std::bitset<3> s)
    {
      for (size_t i=0; i<3; i++)
        sign_[i] = s[i] ? -1.0 : 1.0;
    }

    //! \brief number of shape functions
    unsigned int size () const
    {
      return 6;
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
      out.resize(6);

      out[0][0] = sign_[0]*in[0];
      out[0][1] = sign_[0]*(in[1] - 1.0);
      out[1][0] = sign_[1]*(in[0] - 1.0);
      out[1][1] = sign_[1]*in[1];
      out[2][0] = sign_[2]*in[0];
      out[2][1] = sign_[2]*in[1];
      out[3][0] = 3.0*in[0];
      out[3][1] = 3.0 - 6.0*in[0] - 3.0*in[1];
      out[4][0] = -3.0 + 3.0*in[0] + 6.0*in[1];
      out[4][1] = -3.0*in[1];
      out[5][0] = -3.0*in[0];
      out[5][1] = 3.0*in[1];
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
      out.resize(6);

      out[0][0][0] = sign_[0];
      out[0][0][1] = 0.0;
      out[0][1][0] = 0.0;
      out[0][1][1] = sign_[0];

      out[1][0][0] = sign_[1];
      out[1][0][1] = 0.0;
      out[1][1][0] = 0.0;
      out[1][1][1] = sign_[1];

      out[2][0][0] = sign_[2];
      out[2][0][1] = 0.0;
      out[2][1][0] = 0.0;
      out[2][1][1] = sign_[2];

      out[3][0][0] = 3.0;
      out[3][0][1] = 0.0;
      out[3][1][0] = -6.0;
      out[3][1][1] = -3.0;

      out[4][0][0] = 3.0;
      out[4][0][1] = 6.0;
      out[4][1][0] = 0.0;
      out[4][1][1] = -3.0;

      out[5][0][0] = -3.0;
      out[5][0][1] = 0.0;
      out[5][1][0] = 0.0;
      out[5][1][1] = 3.0;
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
          out[0][0] = sign_[0];
          out[0][1] = 0.0;

          out[1][0] = sign_[1];
          out[1][1] = 0.0;

          out[2][0] = sign_[2];
          out[2][1] = 0.0;

          out[3][0] = 3.0;
          out[3][1] = -6.0;

          out[4][0] = 3.0;
          out[4][1] = 0.0;

          out[5][0] = -3.0;
          out[5][1] = 0.0;
          break;
        case 1:
          out[0][0] = 0.0;
          out[0][1] = sign_[0];

          out[1][0] = 0.0;
          out[1][1] = sign_[1];

          out[2][0] = 0.0;
          out[2][1] = sign_[2];

          out[3][0] = 0.0;
          out[3][1] = -3.0;

          out[4][0] = 6.0;
          out[4][1] = -3.0;

          out[5][0] = 0.0;
          out[5][1] = 3.0;
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
      return 1;
    }

  private:
    std::array<R,3> sign_;
  };
}
#endif // DUNE_LOCALFUNCTIONS_BREZZIDOUGLASMARINI1_SIMPLEX2D_LOCALBASIS_HH
