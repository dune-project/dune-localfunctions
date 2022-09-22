// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_DUAL_Q1_LOCALINTERPOLATION_HH
#define DUNE_DUAL_Q1_LOCALINTERPOLATION_HH

#include <array>
#include <vector>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/localfunctions/common/localinterpolation.hh>

namespace Dune
{

  /** \todo Please doc me! */
  template<int dim, class LB>
  class DualQ1LocalInterpolation
  {
  public:

    void setCoefficients(const std::array<Dune::FieldVector<typename LB::Traits::RangeFieldType, (1<<dim)> ,(1<<dim)>& coefficients)
    {
      coefficients_ = coefficients;
    }


    //! \brief Local interpolation of a function
    template<typename F, typename C>
    void interpolate (const F& ff, std::vector<C>& out) const
    {
      typename LB::Traits::DomainType x;

      auto&& f = Impl::makeFunctionWithCallOperator<decltype(x)>(ff);

      const int size = 1<<dim;

      // compute Q1 interpolation coefficients
      Dune::FieldVector<C,size> q1Coefficients;

      for (int i=0; i< (1<<dim); i++) {

        // Generate coordinate of the i-th corner of the reference cube
        // We could use the ReferenceElement for this as well, but it is
        // still not clear how dune-localfunctions should have access to them.
        for (int j=0; j<dim; j++)
          x[j] = (i & (1<<j)) ? 1.0 : 0.0;

        q1Coefficients[i] = f(x);

      }

      out.resize(size);

      // solve a linear system to compute the dual coefficients
      Dune::FieldMatrix<C,size,size> mat;

      for (int i=0; i<size; i++)
        for (int j=0; j<size; j++)
          mat[i][j] = coefficients_[j][i];

      // now solve for the weights
      Dune::FieldVector<C,size> sol(0);

      mat.solve(sol,q1Coefficients);

      // write result in out vector
      for (int i=0; i<size; i++)
        out[i] = sol[i];
    }

  private:
    std::array<Dune::FieldVector<typename LB::Traits::RangeFieldType, (1<<dim)> ,(1<<dim)> coefficients_;
  };

}

#endif
