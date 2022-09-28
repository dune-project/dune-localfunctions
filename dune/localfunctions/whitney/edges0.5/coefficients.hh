// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

#ifndef DUNE_LOCALFUNCTIONS_WHITNEY_EDGES0_5_COEFFICIENTS_HH
#define DUNE_LOCALFUNCTIONS_WHITNEY_EDGES0_5_COEFFICIENTS_HH

#include <cstddef>
#include <vector>

#include <dune/localfunctions/common/localkey.hh>
#include <dune/localfunctions/whitney/edges0.5/common.hh>

namespace Dune {

  //////////////////////////////////////////////////////////////////////
  //
  // Coefficients
  //

  //! Coefficients for lowest order edge elements on simplices
  /**
   * \nosubgrouping
   * \implements CoefficientsInterface
   *
   * \tparam dim Dimension of both domain and range.
   */
  template<std::size_t dim>
  class EdgeS0_5Coefficients : private EdgeS0_5Common<dim> {
    using EdgeS0_5Common<dim>::s;

    std::vector<LocalKey> li;

  public:
    //! Standard constructor
    EdgeS0_5Coefficients() : li(s) {
      for(std::size_t i = 0; i < s; i++)
        li[i] = LocalKey(i, dim-1, 0);
    }

    //! number of coefficients
    std::size_t size () const { return s; }

    //! get i'th index
    const LocalKey& localKey(std::size_t i) const { return li[i]; }
  };

} // namespace Dune

#endif // DUNE_LOCALFUNCTIONS_WHITNEY_EDGES0_5_COEFFICIENTS_HH
