// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_P0LOCALBASIS_HH
#define DUNE_P0LOCALBASIS_HH

#include <numeric>

#include <dune/common/fmatrix.hh>

#include <dune/localfunctions/common/localbasis.hh>

namespace Dune
{
  /**@ingroup LocalBasisImplementation
         \brief Constant shape function

         Defines the constant scalar shape function in d dimensions. Is
         valid on any type of reference element.

         \tparam D Type to represent the field in the domain.
         \tparam R Type to represent the field in the range.
         \tparam d Domain dimension

         \nosubgrouping
   */
  template<class D, class R, int d>
  class P0LocalBasis
  {
  public:
    //! \brief export type traits for function signature
    typedef LocalBasisTraits<D,d,Dune::FieldVector<D,d>,R,1,Dune::FieldVector<R,1>,
        Dune::FieldMatrix<R,1,d> > Traits;

    //! \brief number of shape functions
    unsigned int size () const
    {
      return 1;
    }

    //! \brief Evaluate all shape functions
    inline void evaluateFunction (const typename Traits::DomainType&,
                                  std::vector<typename Traits::RangeType>& out) const
    {
      out.resize(1);
      out[0] = 1;
    }

    //! \brief Evaluate Jacobian of all shape functions
    inline void
    evaluateJacobian (const typename Traits::DomainType&,         // position
                      std::vector<typename Traits::JacobianType>& out) const      // return value
    {
      out.resize(1);
      for (int i=0; i<d; i++)
        out[0][0][i] = 0;
    }

    /** \brief Evaluate partial derivatives of any order of all shape functions
     * \param order Order of the partial derivatives, in the classic multi-index notation
     * \param in Position where to evaluate the derivatives
     * \param[out] out Return value: the desired partial derivatives
     */
    void partial(const std::array<unsigned int,d>& order,
                 const typename Traits::DomainType& in,
                 std::vector<typename Traits::RangeType>& out) const
    {
      auto totalOrder = std::accumulate(order.begin(), order.end(), 0);
      if (totalOrder == 0) {
        evaluateFunction(in, out);
      } else {
        out.resize(1);
        out[0] = 0;
      }
    }

    //! \brief Polynomial order of the shape functions
    unsigned int order () const
    {
      return 0;
    }
  };

}

#endif
