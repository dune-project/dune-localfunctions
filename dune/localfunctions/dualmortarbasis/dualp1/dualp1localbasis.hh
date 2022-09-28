// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_DUAL_P1_LOCALBASIS_HH
#define DUNE_DUAL_P1_LOCALBASIS_HH

#include <numeric>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/localfunctions/common/localbasis.hh>

namespace Dune
{
  /**@ingroup LocalBasisImplementation
         \brief Dual Lagrange shape functions on the simplex.

         Defines the linear dual shape functions on the simplex.
         Note that if the dual functions are chosen to be dual on the faces,
         the integrated product of a Lagrange \f$\lambda_p\f$ and dual
         function \f$\theta_q\f$ over faces not containing \f$q\f$ does in
         general not vanish.

         \tparam D Type to represent the field in the domain.
         \tparam R Type to represent the field in the range.
         \tparam dim The dimension of the simplex
         \tparam faceDual If set, the basis functions are bi-orthogonal only on faces containing the corresponding vertex.

         \nosubgrouping
   */
  template<class D, class R, int dim, bool faceDualT=false>
  class DualP1LocalBasis
  {
  public:
    //! Determines if the basis is only biorthogonal on adjacent faces
    static const bool faceDual = faceDualT;
    //! \brief export type traits for function signature
    typedef LocalBasisTraits<D,dim,Dune::FieldVector<D,dim>,R,1,Dune::FieldVector<R,1>,
        Dune::FieldMatrix<R,1,dim> > Traits;

    //! \brief number of shape functions
    unsigned int size () const
    {
      return dim+1;
    }

    //! \brief Evaluate all shape functions
    inline void evaluateFunction (const typename Traits::DomainType& in,
                                  std::vector<typename Traits::RangeType>& out) const
    {
      // evaluate P1 basis functions
      std::vector<typename Traits::RangeType> p1Values(size());

      p1Values[0] = 1.0;

      for (int i=0; i<dim; i++) {
        p1Values[0]  -= in[i];
        p1Values[i+1] = in[i];
      }

      // compute dual basis function values as a linear combination of the Lagrange values
      out.resize(size());

      for (int i=0; i<=dim; i++) {
        out[i] = (dim+!faceDual)*p1Values[i];
        for (int j=0; j<i; j++)
          out[i] -= p1Values[j];

        for (int j=i+1; j<=dim; j++)
          out[i] -= p1Values[j];
      }
    }

    //! \brief Evaluate Jacobian of all shape functions
    inline void
    evaluateJacobian (const typename Traits::DomainType& in,
                      std::vector<typename Traits::JacobianType>& out) const
    {
      // evaluate P1 jacobians
      std::vector<typename Traits::JacobianType> p1Jacs(size());

      for (int i=0; i<dim; i++)
        p1Jacs[0][0][i] = -1;

      for (int i=0; i<dim; i++)
        for (int j=0; j<dim; j++)
          p1Jacs[i+1][0][j] = (i==j);

      // compute dual basis jacobians as linear combination of the Lagrange jacobians
      out.resize(size());

      for (size_t i=0; i<=dim; i++) {
        out[i][0] = 0;
        out[i][0].axpy(dim+!faceDual,p1Jacs[i][0]);

        for (size_t j=0; j<i; j++)
          out[i][0] -= p1Jacs[j][0];

        for (int j=i+1; j<=dim; j++)
          out[i][0] -= p1Jacs[j][0];
      }
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
  };
}
#endif
