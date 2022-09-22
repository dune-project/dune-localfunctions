// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_HIERARCHICAL_PRISM_P2_LOCALINTERPOLATION_HH
#define DUNE_HIERARCHICAL_PRISM_P2_LOCALINTERPOLATION_HH

#include <vector>
#include <dune/localfunctions/common/localinterpolation.hh>

namespace Dune
{
  /**
     \tparam LB The LocalBasis implementation
   */
  template<class LB>
  class HierarchicalPrismP2LocalInterpolation
  {
  public:

    template<typename F, typename C>
    void interpolate (const F& ff, std::vector<C>& out) const
    {
      typename LB::Traits::DomainType x;
      typename LB::Traits::RangeType y;
      out.resize(18);

      auto&& f = Impl::makeFunctionWithCallOperator<decltype(x)>(ff);

      //First the  vertex dofs
      x[0] = 0.0;    x[1] = 0.0;     x[2] = 0.0;    out[0] = f(x);
      x[0] = 1.0;    x[1] = 0.0;     x[2] = 0.0;    out[1] = f(x);
      x[0] = 0.0;    x[1] = 1.0;     x[2] = 0.0;    out[2] = f(x);
      x[0] = 0.0;    x[1] = 0.0;     x[2] = 1.0;    out[3] = f(x);
      x[0] = 1.0;    x[1] = 0.0;     x[2] = 1.0;    out[4] = f(x);
      x[0] = 0.0;    x[1] = 1.0;     x[2] = 1.0;    out[5] = f(x);


      // Then: the 9 edge dofs and the 3 face dofs
      x[0] = 0.0;    x[1] = 0.0;     x[2] = 0.5;    y = f(x);
      out[6] = y - 0.5*(out[0] + out[3]);

      x[0] = 1.0;    x[1] = 0.0;     x[2] = 0.5;    y = f(x);
      out[7] = y - 0.5*(out[1] + out[4]);

      x[0] = 0.0;    x[1] = 1.0;     x[2] = 0.5;    y = f(x);
      out[8] = y - 0.5*(out[2] + out[5]);

      x[0] = 0.5;    x[1] = 0.0;     x[2] = 0.0;    y = f(x);
      out[9] = y - 0.5*(out[0] + out[1]);

      x[0] = 0.0;    x[1] = 0.5;     x[2] = 0.0;    y = f(x);
      out[10] = y - 0.5*(out[2] + out[0]);

      x[0] = 0.5;    x[1] = 0.5;     x[2] = 0.0;    y = f(x);
      out[11] = y - 0.5*(out[2] + out[1]);

      x[0] = 0.5;    x[1] = 0.0;     x[2] = 1.0;    y = f(x);
      out[12] = y - 0.5*(out[3] + out[4]);

      x[0] = 0.0;    x[1] = 0.5;     x[2] = 1.0;    y = f(x);
      out[13] = y - 0.5*(out[3] + out[5]);

      x[0] = 0.5;    x[1] = 0.5;     x[2] = 1.0;    y = f(x);
      out[14] = y - 0.5*(out[4] + out[5]);


      //faces
      x[0] = 0.5;    x[1] = 0.0;     x[2] = 0.5;    y = f(x);
      out[15] = y - 0.25*(out[4] + out[1] + out[0] +  out[3] );

      x[0] = 0.0;    x[1] = 0.5;     x[2] = 0.5;    y = f(x);
      out[16] = y - 0.25*(out[2] + out[0] + out[3] + out[5] );

      x[0] = 0.5;    x[1] = 0.5;     x[2] = 0.5;    y = f(x);
      out[17] = y - 0.25*(out[2] + out[1] + out[4] + out[5] );

    }
  };
}

#endif
