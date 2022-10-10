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
   * \ingroup LocalBasisImplementation
   * \brief First order Raviart-Thomas shape functions on the reference pyramid.
   *
   * \tparam D Type to represent the field in the domain.
   * \tparam R Type to represent the field in the range.
   *
   * \nosubgrouping
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
      for (std::size_t i=0; i<out.size(); i++)
        out[i] = {0.0,0.0,0.0};

      out[0][0] = 1.5*in[0];
      out[0][1] = 1.5*in[1];
      out[0][2] = -1.0;

      out[1][0] = -2.0 + 3.0*in[0];

      out[2][0] = 3.0*in[0];

      out[3][1] = -2.0 + 3.0*in[1];

      out[4][1] = 3.0*in[1];

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

      out[0][0][0] = sign[0]*(1.5);
      out[0][1][1] = sign[0]*(1.5);

      out[1][0][0] = sign[1]*(3.0);

      out[2][0][0] = sign[2]*(3.0);

      out[3][1][1] = sign[3]*(3.0);

      out[4][1][1] = sign[4]*(3.0);
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
#endif // DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS0_PYRAMID_LOCALBASIS_HH
