// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS0_PRISM_LOCALBASIS_HH
#define DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS0_PRISM_LOCALBASIS_HH

#include <numeric>
#include <vector>

#include <dune/common/fmatrix.hh>

#include "../../common/localbasis.hh"

namespace Dune
{
  /**
   * \brief First order Raviart-Thomas shape functions on the reference prism.
   *
   * \tparam D Type to represent the field in the domain.
   * \tparam R Type to represent the field in the range.
   *
   * \ingroup RaviartThomasImpl
   * \nosubgrouping
   */
  template<class D, class R>
  class RT0PrismLocalBasis
  {

  public:
    typedef LocalBasisTraits<D,3,Dune::FieldVector<D,3>,R,3,Dune::FieldVector<R,3>,
        Dune::FieldMatrix<R,3,3> > Traits;

    /**
     * \brief Make set number s, where 0 <= s < 32
     *
     * \param s Edge orientation indicator
     */
    RT0PrismLocalBasis (std::bitset<5> s = 0)
    {
      for (size_t i=0; i<size(); i++)
        sign[i] = s[i] ? -1.0 : 1.0;
    }

    //! \brief number of shape functions
    unsigned int size () const
    {
      return 5;
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
      out.resize(5);

      out[0] = {        in[0], -1.0 + in[1],              0.0};

      out[1] = { -1.0 + in[0],        in[1],              0.0};

      out[2] = {        in[0],        in[1],              0.0};

      out[3] = {          0.0,          0.0, -2.0 + 2.0*in[2]};

      out[4] = {          0.0,          0.0,        2.0*in[2]};

      for (std::size_t i=0; i<out.size(); i++)
        out[i] *= sign[i];

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
      out.resize(5);

      for(int i=0; i<size(); i++)
        for(int j=0; j<3; j++)
            out[i][j] = {0.0, 0.0, 0.0};

      out[0][0][0] = sign[0];
      out[0][1][1] = sign[0];

      out[1][0][0] = sign[1];
      out[1][1][1] = sign[1];

      out[2][0][0] = sign[2];
      out[2][1][1] = sign[2];

      out[3][2][2] = sign[3]*(2.0);

      out[4][2][2] = sign[4]*(2.0);
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
    std::array<R,5> sign;
  };
}
#endif // DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS0_PRISM_LOCALBASIS_HH
