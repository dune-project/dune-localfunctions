// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS0_PYRAMID_LOCALINTERPOLATION_HH
#define DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS0_PYRAMID_LOCALINTERPOLATION_HH

#include <vector>

namespace Dune
{
  /**
   * \brief First order Raviart-Thomas shape functions on the reference hexahedron.
   *
   * \tparam LB corresponding LocalBasis giving traits
   *
   * \nosubgrouping
   * \ingroup RaviartThomasImpl
   */
  template<class LB>
  class RT0PyramidLocalInterpolation
  {

  public:

    /**
     * \brief Make set number s, where 0 <= s < 32
     *
     * \param s Face orientation indicator
     */
    RT0PyramidLocalInterpolation (std::bitset<5> s = 0)
    {
      typedef typename LB::Traits::RangeFieldType Scalar;

      for (size_t i=0; i<5; i++)
        sign_[i] = (s[i]) ? -1.0 : 1.0;

      // No need to flip the sign for the interior basis function
      sign_[5] = 1.0;

      Scalar r = 1/std::sqrt(2);

      facetNormal_[0] = { 0.0,  0.0, -1.0};
      facetNormal_[1] = {-1.0,  0.0,  0.0};
      facetNormal_[2] = {   r,  0.0,    r};
      facetNormal_[3] = { 0.0, -1.0,  0.0};
      facetNormal_[4] = { 0.0,    r,    r};
      facetNormal_[5] = {   r,   -r,  0.0};

      facetArea_[0] = 1.0;
      facetArea_[1] = 1/2.0;
      facetArea_[2] = 1/2.0 * std::sqrt(2);
      facetArea_[3] = 1/2.0;
      facetArea_[4] = 1/2.0 * std::sqrt(2);
      facetArea_[5] = 1/2.0 * std::sqrt(2);

      facetCenter_[0] = {   0.5,   0.5,   0.0};
      facetCenter_[1] = {   0.0, 1/3.0, 1/3.0};
      facetCenter_[2] = { 2/3.0, 1/3.0, 1/3.0};
      facetCenter_[3] = { 1/3.0,   0.0, 1/3.0};
      facetCenter_[4] = { 1/3.0, 2/3.0, 1/3.0};
      facetCenter_[5] = { 1/3.0, 1/3.0, 1/3.0};
    }

    /**
     * \brief Interpolate a given function with shape functions
     *
     * \tparam F Function type for function which should be interpolated
     * \tparam C Coefficient type
     * \param f function which should be interpolated
     * \param out return value, vector of coefficients
     */
    template<class F, class C>
    void interpolate (const F& f, std::vector<C>& out) const
    {
      out.resize(6);
      for(int i=0; i<6; i++)
        out[i] = f(facetCenter_[i]).dot(facetNormal_[i]) * facetArea_[i] * sign_[i];

      // Adjust the scaling for the triangular faces functions.
      // The tetrahedral RT-basis is scaled rather strange:
      // The dual basis does not evaluate the face integral
      // of the normal component, but sqrt(2) times this integral.
      // In order to match the basis functions, we need to rescale
      // the triangular face functions here, too.
      for(std::size_t i=1; i<5; ++i)
        out[i] *= std::sqrt(2.0);
    }

  private:
    // Facet orientations
    std::array<typename LB::Traits::RangeFieldType, 6> sign_;
    // Facet area
    std::array<typename LB::Traits::RangeFieldType, 6> facetArea_;

    // Facet normals
    std::array<typename LB::Traits::DomainType, 6> facetNormal_;
    // Facet midpoints
    std::array<typename LB::Traits::DomainType, 6> facetCenter_;
  };
}
#endif // DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS0_PYRAMID_LOCALINTERPOLATION_HH
