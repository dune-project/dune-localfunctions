// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_HIERARCHICAL_SIMPLEX_P2_LOCALINTERPOLATION_HH
#define DUNE_HIERARCHICAL_SIMPLEX_P2_LOCALINTERPOLATION_HH

#include <vector>
#include <dune/localfunctions/common/localinterpolation.hh>

namespace Dune
{
  /**
     \tparam LB The LocalBasis implementation
   */
  template<class LB>
  class HierarchicalSimplexP2LocalInterpolation
  {
  public:

    template<typename F, typename C>
    void interpolate (const F& ff, std::vector<C>& out) const
    {
      typename LB::Traits::DomainType x;
      typename LB::Traits::RangeType y;

      auto&& f = Impl::makeFunctionWithCallOperator<decltype(x)>(ff);

      static_assert(LB::Traits::dimDomain <= 3,
                    "LocalInterpolation for HierarchicalSimplexP2 finite elements"
                    " is only implemented for dimDomain <=3!");

      switch ( int(LB::Traits::dimDomain))  {

      case 1 :

        out.resize(3);

        // First: the two vertex dofs
        x[0] = 0.0;   out[0] = f(x);
        x[0] = 1.0;   out[2] = f(x);

        // Then: the edge dof
        x[0] = 0.5;   y = f(x);
        out[1] = y - 0.5*(out[0] + out[2]);

        break;


      case 2 :

        out.resize(6);

        // First: the three vertex dofs
        x[0] = 0.0;    x[1] = 0.0;      out[0] = f(x);
        x[0] = 1.0;    x[1] = 0.0;      out[2] = f(x);
        x[0] = 0.0;    x[1] = 1.0;      out[5] = f(x);

        // Then: the three edge dofs
        x[0] = 0.5;    x[1] = 0.0;      y = f(x);
        out[1] = y - 0.5*(out[0] + out[2]);

        x[0] = 0.0;    x[1] = 0.5;      y = f(x);
        out[3] = y - 0.5*(out[0] + out[5]);

        x[0] = 0.5;    x[1] = 0.5;      y = f(x);
        out[4] = y - 0.5*(out[2] + out[5]);

        break;

      case 3 :

        out.resize(10);

        // First: the four vertex dofs
        x[0] = 0.0;    x[1] = 0.0;     x[2] = 0.0;    out[0] = f(x);
        x[0] = 1.0;    x[1] = 0.0;     x[2] = 0.0;    out[2] = f(x);
        x[0] = 0.0;    x[1] = 1.0;     x[2] = 0.0;    out[5] = f(x);
        x[0] = 0.0;    x[1] = 0.0;     x[2] = 1.0;    out[9] = f(x);

        // Then: the six edge dofs
        x[0] = 0.5;    x[1] = 0.0;     x[2] = 0.0;    y = f(x);
        out[1] = y - 0.5*(out[0] + out[2]);

        x[0] = 0.0;    x[1] = 0.5;     x[2] = 0.0;    y = f(x);
        out[3] = y - 0.5*(out[0] + out[5]);

        x[0] = 0.5;    x[1] = 0.5;     x[2] = 0.0;    y = f(x);
        out[4] = y - 0.5*(out[2] + out[5]);

        x[0] = 0.0;    x[1] = 0.0;     x[2] = 0.5;    y = f(x);
        out[6] = y - 0.5*(out[0] + out[9]);

        x[0] = 0.5;    x[1] = 0.0;     x[2] = 0.5;    y = f(x);
        out[7] = y - 0.5*(out[2] + out[9]);

        x[0] = 0.0;    x[1] = 0.5;     x[2] = 0.5;    y = f(x);
        out[8] = y - 0.5*(out[5] + out[9]);

        break;

      }
    }

  };
}

#endif
