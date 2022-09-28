// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_DUAL_P1_LOCALINTERPOLATION_HH
#define DUNE_DUAL_P1_LOCALINTERPOLATION_HH

#include <vector>
#include <dune/localfunctions/common/localinterpolation.hh>

namespace Dune
{
  template<int dim, class LB>
  class DualP1LocalInterpolation
  {
  public:
    //! \brief Local interpolation of a function
    template<typename F, typename C>
    void interpolate (const F& ff, std::vector<C>& out) const
    {
      typename LB::Traits::DomainType x;
      // If the dual functions are dual on the faces,
      // then adjust the interpolation weights
      const int faceDual(LB::faceDual);

      auto&& f = Impl::makeFunctionWithCallOperator<decltype(x)>(ff);

      // compute P1 interpolation coefficients
      std::vector<C> p1Interpolation(dim+1);

      // vertex 0
      for (int i=0; i<dim; i++)
        x[i] = 0;
      p1Interpolation[0] = f(x);

      // remaining vertices
      for (int i=0; i<dim; i++) {
        for (int j=0; j<dim; j++)
          x[j] = (i==j);

        p1Interpolation[i+1] = f(x);

      }

      // compute dual coefficients from the Lagrange ones
      out.resize(dim+1);
      for (int i=0; i<dim+1; i++) {
        out[i] = 2*p1Interpolation[i]/(dim+2-faceDual);

        for (int j=0; j<i; j++)
          out[i] += p1Interpolation[j]/(dim+2-faceDual);

        for (int j=i+1; j<=dim; j++)
          out[i] += p1Interpolation[j]/(dim+2-faceDual);
      }
    }

  };
}

#endif
