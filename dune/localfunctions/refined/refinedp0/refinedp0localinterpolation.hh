// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_REFINED_P0_LOCALINTERPOLATION_HH
#define DUNE_REFINED_P0_LOCALINTERPOLATION_HH

#include <dune/localfunctions/refined/refinedp0/refinedp0localbasis.hh>
#include <dune/localfunctions/common/localinterpolation.hh>

namespace Dune
{
  template<class LB>
  class RefinedP0LocalInterpolation
  {};

  template<class D, class R>
  class RefinedP0LocalInterpolation<RefinedP0LocalBasis<D,R,1> >
  {
    typedef RefinedP0LocalBasis<D,R,1> LB;
    typedef typename LB::Traits::DomainType DT;

  public:
    RefinedP0LocalInterpolation() :
      interpolationPoints_(2)
    {
      // Interpolation is done by evaluating at the halved segments centers
      interpolationPoints_[0][0] = 1.0/4.0;

      interpolationPoints_[1][0] = 3.0/4.0;
    }


    template<typename F, typename C>
    void interpolate (const F& ff, std::vector<C>& out) const
    {
      auto&& f = Impl::makeFunctionWithCallOperator<DT>(ff);

      out.resize(interpolationPoints_.size());
      for (size_t i = 0; i < out.size(); ++i)
      {
        out[i] = f(interpolationPoints_[i]);
      }
    }

  private:
    std::vector<DT> interpolationPoints_;
  };

  template<class D, class R>
  class RefinedP0LocalInterpolation<RefinedP0LocalBasis<D,R,2> >
  {
    typedef RefinedP0LocalBasis<D,R,2> LB;
    typedef typename LB::Traits::DomainType DT;

  public:
    RefinedP0LocalInterpolation() :
      interpolationPoints_(4)
    {
      // Interpolation is done by evaluating at the subtriangle centers
      interpolationPoints_[0][0] = 1.0/6;
      interpolationPoints_[0][1] = 1.0/6;

      interpolationPoints_[1][0] = 4.0/6;
      interpolationPoints_[1][1] = 1.0/6;

      interpolationPoints_[2][0] = 1.0/6;
      interpolationPoints_[2][1] = 4.0/6;

      interpolationPoints_[3][0] = 2.0/6;
      interpolationPoints_[3][1] = 2.0/6;
    }


    template<typename F, typename C>
    void interpolate (const F& ff, std::vector<C>& out) const
    {
      auto&& f = Impl::makeFunctionWithCallOperator<DT>(ff);

      out.resize(interpolationPoints_.size());
      for (size_t i = 0; i < out.size(); ++i)
      {
        out[i] = f(interpolationPoints_[i]);
      }
    }

  private:
    std::vector<DT> interpolationPoints_;
  };

  template<class D, class R>
  class RefinedP0LocalInterpolation<RefinedP0LocalBasis<D,R,3> >
  {
    typedef RefinedP0LocalBasis<D,R,3> LB;
    typedef typename LB::Traits::DomainType DT;

  public:
    RefinedP0LocalInterpolation() :
      interpolationPoints_(8)
    {
      // Interpolation is done by evaluating at the subtriangle centers
      interpolationPoints_[0][0] = 1.0/8;
      interpolationPoints_[0][1] = 1.0/8;
      interpolationPoints_[0][2] = 1.0/8;

      interpolationPoints_[1][0] = 5.0/8;
      interpolationPoints_[1][1] = 1.0/8;
      interpolationPoints_[1][2] = 1.0/8;

      interpolationPoints_[2][0] = 1.0/8;
      interpolationPoints_[2][1] = 5.0/8;
      interpolationPoints_[2][2] = 1.0/8;

      interpolationPoints_[3][0] = 1.0/8;
      interpolationPoints_[3][1] = 1.0/8;
      interpolationPoints_[3][2] = 5.0/8;

      interpolationPoints_[4][0] = 1.0/4;
      interpolationPoints_[4][1] = 1.0/8;
      interpolationPoints_[4][2] = 1.0/4;

      interpolationPoints_[5][0] = 3.0/8;
      interpolationPoints_[5][1] = 1.0/4;
      interpolationPoints_[5][2] = 1.0/8;

      interpolationPoints_[6][0] = 1.0/8;
      interpolationPoints_[6][1] = 1.0/4;
      interpolationPoints_[6][2] = 3.0/8;

      interpolationPoints_[7][0] = 1.0/4;
      interpolationPoints_[7][1] = 3.0/8;
      interpolationPoints_[7][2] = 1.0/4;
    }


    template<typename F, typename C>
    void interpolate (const F& ff, std::vector<C>& out) const
    {
      auto&& f = Impl::makeFunctionWithCallOperator<DT>(ff);

      out.resize(interpolationPoints_.size());
      for (size_t i = 0; i < out.size(); ++i)
      {
        out[i] = f(interpolationPoints_[i]);
      }
    }

  private:
    std::vector<DT> interpolationPoints_;
  };
}

#endif
