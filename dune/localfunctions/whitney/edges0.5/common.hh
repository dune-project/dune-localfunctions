// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_LOCALFUNCTIONS_WHITNEY_EDGES0_5_COMMON_HH
#define DUNE_LOCALFUNCTIONS_WHITNEY_EDGES0_5_COMMON_HH

#include <cstddef>

#include <dune/geometry/referenceelements.hh>

namespace Dune {

  //! Common base class for edge elements
  template<std::size_t dim, class DF = double>
  struct EdgeS0_5Common {
    //! The reference element for this edge element
    static const ReferenceElement<DF, dim>& refelem;
    //! The number of base functions
    /**
     * \note This is not a compile time constant, since the number of edges is
     *       extracted from the reference element.
     */
    static const std::size_t s;
  };

  template<std::size_t dim, class DF>
  const ReferenceElement<DF, dim>& EdgeS0_5Common<dim,DF>::
  refelem(ReferenceElements<DF, dim>::simplex());

  template<std::size_t dim, typename DF>
  const std::size_t EdgeS0_5Common<dim,DF>::s(refelem.size(dim-1));

} // namespace Dune

#endif // DUNE_LOCALFUNCTIONS_WHITNEY_EDGES0_5_COMMON_HH
