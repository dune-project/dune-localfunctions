// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_ENRICHED_SIMPLEXP1BUBBLE_LOCALCOEFFICIENTS_HH
#define DUNE_LOCALFUNCTIONS_ENRICHED_SIMPLEXP1BUBBLE_LOCALCOEFFICIENTS_HH

#include <array>

#include <dune/localfunctions/common/localkey.hh>

namespace Dune
{
  /**
   * \ingroup LocalBasisImplementation
   * \brief The Local keys associated to the dim-d local basis functions
   *
   * The local keys represent first the vertex DOFs and then a single
   * DOF associated to the element interior.
   * \nosubgrouping
   **/
  template <int dim>
  class SimplexP1BubbleLocalCoefficients
  {
    static const int numVertices = dim+1;

  public:
    //! Default constructor, initialized the local keys
    SimplexP1BubbleLocalCoefficients () noexcept
    {
      for (int i = 0; i <= dim; ++i)
        li_[i] = LocalKey(i,dim,0);       // Vertex
      li_[numVertices] = LocalKey(0,0,0); // Element
    }

    //! Returns number of coefficients
    static constexpr std::size_t size () noexcept
    {
      return numVertices + 1;
    }

    //! Returns the i'th local key
    const LocalKey& localKey (std::size_t i) const noexcept
    {
      return li_[i];
    }

  private:
    std::array<LocalKey, numVertices+1> li_;
  };

} // end namespace Dune

#endif // DUNE_LOCALFUNCTIONS_ENRICHED_SIMPLEXP1BUBBLE_LOCALCOEFFICIENTS_HH
