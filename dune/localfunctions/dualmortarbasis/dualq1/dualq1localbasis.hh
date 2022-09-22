// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_DUAL_Q1_LOCALBASIS_HH
#define DUNE_DUAL_Q1_LOCALBASIS_HH

#include <array>
#include <numeric>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <dune/localfunctions/common/localbasis.hh>

namespace Dune
{
  /**@ingroup LocalBasisImplementation
         \brief Dual Lagrange shape functions of order 1 on the reference cube.

         \tparam D Type to represent the field in the domain.
         \tparam R Type to represent the field in the range.
     \tparam dim Dimension of the cube

         \nosubgrouping
   */
  template<class D, class R, int dim>
  class DualQ1LocalBasis
  {
  public:
    typedef LocalBasisTraits<D,dim,Dune::FieldVector<D,dim>,R,1,Dune::FieldVector<R,1>,
        Dune::FieldMatrix<R,1,dim> > Traits;

    void setCoefficients(const std::array<Dune::FieldVector<R, (1<<dim)> ,(1<<dim)>& coefficients)
    {
      coefficients_ = coefficients;
    }

    //! \brief number of shape functions
    unsigned int size () const
    {
      return 1<<dim;
    }

    //! \brief Evaluate all shape functions
    inline void evaluateFunction (const typename Traits::DomainType& in,
                                  std::vector<typename Traits::RangeType>& out) const
    {
      // compute q1 values
      std::vector<typename Traits::RangeType> q1Values(size());

      for (size_t i=0; i<size(); i++) {

        q1Values[i] = 1;

        for (int j=0; j<dim; j++)
          // if j-th bit of i is set multiply with in[j], else with 1-in[j]
          q1Values[i] *= (i & (1<<j)) ? in[j] :  1-in[j];

      }

      // compute the dual values by using that they are linear combinations of q1 functions
      out.resize(size());
      for (size_t i=0; i<size(); i++)
        out[i] = 0;

      for (size_t i=0; i<size(); i++)
        for (size_t j=0; j<size(); j++)
          out[i] += coefficients_[i][j]*q1Values[j];


    }

    //! \brief Evaluate Jacobian of all shape functions
    inline void
    evaluateJacobian (const typename Traits::DomainType& in,             // position
                      std::vector<typename Traits::JacobianType>& out) const // return value
    {
      // compute q1 jacobians
      std::vector<typename Traits::JacobianType> q1Jacs(size());

      // Loop over all shape functions
      for (size_t i=0; i<size(); i++) {

        // Loop over all coordinate directions
        for (int j=0; j<dim; j++) {

          // Initialize: the overall expression is a product
          // if j-th bit of i is set to -1, else 1
          q1Jacs[i][0][j] = (i & (1<<j)) ? 1 : -1;

          for (int k=0; k<dim; k++) {

            if (j!=k)
              // if k-th bit of i is set multiply with in[j], else with 1-in[j]
              q1Jacs[i][0][j] *= (i & (1<<k)) ? in[k] :  1-in[k];

          }

        }

      }

      // compute the dual jacobians by using that they are linear combinations of q1 functions
      out.resize(size());
      for (size_t i=0; i<size(); i++)
        out[i] = 0;

      for (size_t i=0; i<size(); i++)
        for (size_t j=0; j<size(); j++)
          out[i].axpy(coefficients_[i][j],q1Jacs[j]);

    }

    //! \brief Evaluate partial derivatives of all shape functions
    void partial (const std::array<unsigned int, dim>& order,
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
    std::array<Dune::FieldVector<R, (1<<dim)> ,(1<<dim)> coefficients_;
  };
}
#endif
