// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_P0LOCALINTERPOLATION_HH
#define DUNE_P0LOCALINTERPOLATION_HH

#include <vector>
#include <dune/geometry/referenceelements.hh>
#include <dune/localfunctions/common/localinterpolation.hh>


namespace Dune
{

  template<class LB>
  class P0LocalInterpolation
  {
  public:
    P0LocalInterpolation (const GeometryType& gt) : gt_(gt)
    {}

    //! determine coefficients interpolating a given function
    template<typename F, typename C>
    void interpolate (const F& ff, std::vector<C>& out) const
    {
      typedef typename LB::Traits::DomainType DomainType;
      typedef typename LB::Traits::DomainFieldType DF;
      const int dim=LB::Traits::dimDomain;

      auto&& f = Impl::makeFunctionWithCallOperator<typename LB::Traits::DomainType>(ff);

      DomainType x = Dune::ReferenceElements<DF,dim>::general(gt_).position(0,0);

      out.resize(1);
      out[0] = f(x);
    }

  private:
    GeometryType gt_;
  };

}

#endif
