// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

#ifndef DUNE_LOCALFUNCTIONS_WHITNEY_EDGES0_5_COMMON_HH
#define DUNE_LOCALFUNCTIONS_WHITNEY_EDGES0_5_COMMON_HH

#include <cstddef>

#include <dune/geometry/dimension.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>

namespace Dune {

  //! Common base class for edge elements
  template<std::size_t dim, class DF = double>
  struct EdgeS0_5Common {
    //! The type of the referenceElement
    using RefElem =
      decltype(referenceElement(DF{}, GeometryTypes::simplex(dim),
                                Dim<dim>{}));

    //! The reference element for this edge element
    RefElem refelem = referenceElement(DF{}, GeometryTypes::simplex(dim),
                                       Dim<dim>{});

    //! The number of base functions
    /**
     * \note This is not a compile time constant, since the number of edges is
     *       extracted from the reference element.
     */
    std::size_t s = refelem.size(dim-1);
  };

} // namespace Dune

#endif // DUNE_LOCALFUNCTIONS_WHITNEY_EDGES0_5_COMMON_HH
