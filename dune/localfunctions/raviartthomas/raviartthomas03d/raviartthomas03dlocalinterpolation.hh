// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS_RAVIARTTHOMAS03D_RAVIARTTHOMAS03DLOCALINTERPOLATION_HH
#define DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS_RAVIARTTHOMAS03D_RAVIARTTHOMAS03DLOCALINTERPOLATION_HH

#include <cmath>
#include <array>
#include <bitset>
#include <vector>
#include <dune/localfunctions/common/localinterpolation.hh>

namespace Dune
{
  template<class LB>
  class RT03DLocalInterpolation
  {
  public:

    //! \brief Constructor with given set of face orientations
    RT03DLocalInterpolation (std::bitset<4> s = 0)
    {
      using std::sqrt;
      for (std::size_t i=0; i<sign_.size(); i++)
        sign_[i] = (s[i]) ? -1.0 : 1.0;

      m_[0] = {1/3.0, 1/3.0,   0.0};
      m_[1] = {1/3.0,   0.0, 1/3.0};
      m_[2] = {  0.0, 1/3.0, 1/3.0};
      m_[3] = {1/3.0, 1/3.0, 1/3.0};
      n_[0] = {          0.0,           0.0,          -1.0};
      n_[1] = {          0.0,          -1.0,           0.0};
      n_[2] = {         -1.0,           0.0,           0.0};
      n_[3] = {1.0/sqrt(3.0), 1.0/sqrt(3.0), 1.0/sqrt(3.0)};
      c_[0] = sqrt(2.0);
      c_[1] = sqrt(2.0);
      c_[2] = sqrt(2.0);
      c_[3] = sqrt(2.0)/sqrt(3.0);
    }

    template<typename F, typename C>
    void interpolate (const F& ff, std::vector<C>& out) const
    {
      // f gives v*outer normal at a point on the face!
      auto&& f = Impl::makeFunctionWithCallOperator<typename LB::Traits::DomainType>(ff);

      out.resize(4);

      for (int i=0; i<4; i++)
      {
        auto y = f(m_[i]);
        out[i] = (y[0]*n_[i][0]+y[1]*n_[i][1]+y[2]*n_[i][2])*sign_[i]/c_[i];
      }
    }

  private:
    // Face orientations
    std::array<typename LB::Traits::RangeFieldType,4> sign_;
    // Face midpoints of the reference tetrahedron
    std::array<typename LB::Traits::DomainType,4> m_;
    // Unit outer normals of the reference tetrahedron
    std::array<typename LB::Traits::DomainType,4> n_;
    // Inverse triangle face area
    std::array<typename LB::Traits::RangeFieldType,4> c_;
  };
}

#endif
