// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS0_PYRAMID_LOCALBASIS_HH
#define DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS0_PYRAMID_LOCALBASIS_HH

#include <numeric>
#include <vector>

#include <dune/common/fmatrix.hh>

#include "../../common/localbasis.hh"

namespace Dune
{
  /**
   * \brief First order Raviart-Thomas shape functions on the reference pyramid.
   *
   * \tparam D Type to represent the field in the domain.
   * \tparam R Type to represent the field in the range.
   *
   * \nosubgrouping
   * \ingroup RaviartThomasImpl
   */
  template<class D, class R>
  class RT0PyramidLocalBasis
  {

  public:
    typedef LocalBasisTraits<D,3,Dune::FieldVector<D,3>,R,3,Dune::FieldVector<R,3>,
        Dune::FieldMatrix<R,3,3> > Traits;

    /**
     * \brief Make set number s, where 0 <= s < 32
     *
     * \param s Edge orientation indicator
     */
    RT0PyramidLocalBasis (std::bitset<5> s = 0)
    {
      // For each basis function we store the constant offset vector
      // and the factor of the linear term for both sub-elements.
      // In the j-th sub-element the i-th basis function is then
      // given by y = offset_[i][j] + factor_[i][j]*x.
      offset_[0][0] = {0.0, 0.0, -1.0};
      factor_[0][0] = 1.0;
      offset_[0][1] = {0.0, 0.0, -1.0};
      factor_[0][1] = 1.0;

      offset_[1][0] = {-2.0, -2.0, 0.0};
      factor_[1][0] = 2.0;
      offset_[1][1] = {0.0, 0.0, 0.0};
      factor_[1][1] = 0.0;

      offset_[2][0] = {0.0, 0.0, 0.0};
      factor_[2][0] = 0.0;
      offset_[2][1] = {0.0, 0.0, 0.0};
      factor_[2][1] = 2.0;

      offset_[3][0] = {0.0, 0.0, 0.0};
      factor_[3][0] = 0.0;
      offset_[3][1] = {-2.0, -2.0, 0.0};
      factor_[3][1] = 2.0;

      offset_[4][0] = {0.0, 0.0, 0.0};
      factor_[4][0] = 2.0;
      offset_[4][1] = {0.0, 0.0, 0.0};
      factor_[4][1] = 0.0;

      // Interior basis function associated to the interior
      // face given by the intersection of the sub-elements
      // in the plane where x[0]==x[1].
      offset_[5][0] = {0.0, -2.0, 0.0};
      factor_[5][0] = 2.0;
      offset_[5][1] = {2.0, 0.0, 0.0};
      factor_[5][1] = -2.0;

      for (size_t i=0; i<5; i++)
        sign_[i] = s[i] ? -1.0 : 1.0;

      // No need to flip the sign_ for the interior basis function
      sign_[5] = 1.0;
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
    inline void evaluateFunction (const typename Traits::DomainType& x,
                                  std::vector<typename Traits::RangeType>& out) const
    {
      out.resize(size());
      bool compositeElement = x[0] > x[1];
      for (std::size_t i=0; i<size(); i++)
      {
        out[i] = x;
        out[i] *= factor_[i][compositeElement];
        out[i] += offset_[i][compositeElement];
        out[i] *= sign_[i];
      }
    }

    /**
     * \brief Evaluate Jacobian of all shape functions
     *
     * \param in Position
     * \param out return value
     */
    inline void evaluateJacobian (const typename Traits::DomainType& x,
                                  std::vector<typename Traits::JacobianType>& out) const
    {
      out.resize(size());
      bool compositeElement = x[0] > x[1];
      for (std::size_t i=0; i<size(); i++)
        for(std::size_t j=0; j<3; j++)
          for(std::size_t k=0; k<3; k++)
            out[i][j][k] = (j==k) * factor_[i][compositeElement] * sign_[i];
    }

    //! \brief Evaluate partial derivatives of all shape functions
    void partial (const std::array<unsigned int, 3>& order,
                  const typename Traits::DomainType& in,         // position
                  std::vector<typename Traits::RangeType>& out) const      // return value
    {
      auto totalOrder = std::accumulate(order.begin(), order.end(), 0);
      if (totalOrder == 0) {
        evaluateFunction(in, out);
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
    std::array<R,6> sign_;
    std::array<std::array<typename Traits::RangeType, 2>, 6> offset_;
    std::array<std::array<R, 2>, 6> factor_;
  };
}
#endif // DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS0_PYRAMID_LOCALBASIS_HH
