// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LOCALFUNCTIONS_COMMON_PARTIAL_HH
#define DUNE_LOCALFUNCTIONS_COMMON_PARTIAL_HH

#include <array>
#include <cstddef>

namespace Dune
{
  namespace Impl
  {
    // helper function to convert an differentiation-order array used in the
    // LocalBasis::partial methods to a differentiation-direction array used in
    // the LocalBasis::evaluate methods.
    template <std::size_t dim, std::size_t diffOrder>
    void order2directions(std::array<unsigned int, dim> const& order, std::array<int, diffOrder>& directions)
    {
      std::size_t counter = 0;
      for (int i = 0; i < int(dim); ++i)
        for (unsigned int j = 0; j < order[i]; ++j)
          directions[counter++] = i;
    }

    // helper function to convert an differentiation-direction array used in the
    // LocalBasis::evaluate methods to a differentiation-order array used in
    // the LocalBasis::partial methods.
    template <std::size_t dim, std::size_t diffOrder>
    void directions2order(std::array<int, diffOrder> const& directions, std::array<unsigned int, dim>& order)
    {
      for(std::size_t i = 0; i < dim; ++i)
        order[i] = 0;
      for(std::size_t i = 0; i < diffOrder; ++i)
        order[directions[i]]++;
    }

  } // end namespace Impl
} // namespace Dune

#endif //DUNE_LOCALFUNCTIONS_COMMON_PARTIAL_HH
