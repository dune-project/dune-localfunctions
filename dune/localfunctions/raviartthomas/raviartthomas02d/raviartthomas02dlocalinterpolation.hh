// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_RT02DLOCALINTERPOLATION_HH
#define DUNE_RT02DLOCALINTERPOLATION_HH

#include <cmath>
#include <array>
#include <bitset>
#include <vector>
#include <dune/localfunctions/common/localinterpolation.hh>

namespace Dune
{
  template<class LB>
  class RT02DLocalInterpolation
  {
  public:

    //! \brief Constructor with given set of edge orientations
    RT02DLocalInterpolation (std::bitset<3> s = 0)
    {
      using std::sqrt;
      for (std::size_t i=0; i<sign_.size(); i++)
        sign_[i] = (s[i]) ? -1.0 : 1.0;

      m_[0] = {0.5, 0.0};
      m_[1] = {0.0, 0.5};
      m_[2] = {0.5, 0.5};
      n_[0] = {0.0,          -1.0};
      n_[1] = {-1.0,          0.0};
      n_[2] = {1.0/sqrt(2.0), 1.0/sqrt(2.0)};
      c_[0] = ( 0.5*n_[0][0] - 1.0*n_[0][1]);
      c_[1] = (-1.0*n_[1][0] + 0.5*n_[1][1]);
      c_[2] = ( 0.5*n_[2][0] + 0.5*n_[2][1]);
    }

    template<typename F, typename C>
    void interpolate (const F& ff, std::vector<C>& out) const
    {
      // f gives v*outer normal at a point on the edge!
      auto&& f = Impl::makeFunctionWithCallOperator<typename LB::Traits::DomainType>(ff);

      out.resize(3);

      for (int i=0; i<3; i++)
      {
        auto y = f(m_[i]);
        out[i] = (y[0]*n_[i][0]+y[1]*n_[i][1])*sign_[i]/c_[i];
      }
    }

  private:
    // Edge orientations
    std::array<typename LB::Traits::RangeFieldType,3> sign_;
    // Edge midpoints of the reference triangle
    std::array<typename LB::Traits::DomainType,3> m_;
    // Unit outer normals of the reference triangle
    std::array<typename LB::Traits::DomainType,3> n_;
    // Inverse triangle edge length
    std::array<typename LB::Traits::RangeFieldType,3> c_;
  };
}

#endif
