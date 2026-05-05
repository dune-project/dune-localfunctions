// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_LAGRANGE_LAGRANGELFECACHE_HH
#define DUNE_LOCALFUNCTIONS_LAGRANGE_LAGRANGELFECACHE_HH

#include <tuple>
#include <utility>

#include <dune/geometry/type.hh>
#include <dune/geometry/typeindex.hh>

#include <dune/localfunctions/lagrange/lagrangecube.hh>
#include <dune/localfunctions/lagrange/lagrangeprism.hh>
#include <dune/localfunctions/lagrange/lagrangepyramid.hh>
#include <dune/localfunctions/lagrange/lagrangesimplex.hh>
#include <dune/localfunctions/lagrange/p0.hh>
#include <dune/localfunctions/common/localfiniteelementvariantcache.hh>


namespace Dune {



namespace Impl {

  // Provide implemented Lagrange local finite elements

  template<class D, class R, std::size_t dim, int compileTimeOrder>
  class ImplementedLagrangeFiniteElements : public FixedDimLocalGeometryTypeIndex<dim>
  {
    int runTimeOrder_;
  public:

    ImplementedLagrangeFiniteElements(int runTimeOrder = compileTimeOrder)
      : runTimeOrder_(runTimeOrder)
    {
      if (runTimeOrder < 0)
        DUNE_THROW(Dune::InvalidStateException, "LagrangeLocalFiniteElementCache: Run-time order must be non-negative.");
      if ((compileTimeOrder >= 0) and (runTimeOrder != compileTimeOrder))
        DUNE_THROW(Dune::InvalidStateException, "LagrangeLocalFiniteElementCache: Run-time order must be consistent with compile-time order when providing both.");
    }

    using FixedDimLocalGeometryTypeIndex<dim>::index;

    auto getImplementations() const
    {
      if constexpr (compileTimeOrder >= 0)
      {
        if constexpr ((dim != 3) and (compileTimeOrder == 0))
          return std::make_tuple(
            std::make_pair(index(GeometryTypes::simplex(dim)), []() { return P0LocalFiniteElement<D,R,dim>(GeometryTypes::simplex(dim)); }),
            std::make_pair(index(GeometryTypes::cube(dim)),    []() { return P0LocalFiniteElement<D,R,dim>(GeometryTypes::cube(dim)); }),
            std::make_pair(index(GeometryTypes::none(dim)),    []() { return P0LocalFiniteElement<D,R,dim>(GeometryTypes::none(dim)); })
          );
        else if constexpr ((dim == 3) and (compileTimeOrder <= 2))
          return std::make_tuple(
            std::make_pair(index(GeometryTypes::tetrahedron), []() { return LagrangeSimplexLocalFiniteElement<D,R,dim,compileTimeOrder>(); }),
            std::make_pair(index(GeometryTypes::hexahedron),  []() { return LagrangeCubeLocalFiniteElement<D,R,dim,compileTimeOrder>(); }),
            std::make_pair(index(GeometryTypes::prism),       []() { return LagrangePrismLocalFiniteElement<D,R,compileTimeOrder>(); }),
            std::make_pair(index(GeometryTypes::pyramid),     []() { return LagrangePyramidLocalFiniteElement<D,R,compileTimeOrder>(); })
          );
        else
          return std::make_tuple(
            std::make_pair(index(GeometryTypes::simplex(dim)), []() { return LagrangeSimplexLocalFiniteElement<D,R,dim,compileTimeOrder>(); }),
            std::make_pair(index(GeometryTypes::cube(dim)),    []() { return LagrangeCubeLocalFiniteElement<D,R,dim,compileTimeOrder>(); })
          );
      }
      else
      {
        if constexpr (dim == 3)
        {
          constexpr auto unusedIndex = std::numeric_limits<decltype(index(GeometryTypes::simplex(dim)))>::max();
          auto prismOrder = (runTimeOrder_<= 2) ? runTimeOrder_ : 0;
          auto prismIndex = (runTimeOrder_<= 2) ? index(GeometryTypes::prism) : unusedIndex;
          auto pyramidOrder = (runTimeOrder_<= 2) ? runTimeOrder_ : 0;
          auto pyramidIndex = (runTimeOrder_<= 2) ? index(GeometryTypes::pyramid) : unusedIndex;
          return std::make_tuple(
            std::make_pair(index(GeometryTypes::tetrahedron), [&]() { return LagrangeSimplexLocalFiniteElement<D,R,dim,compileTimeOrder>(runTimeOrder_); }),
            std::make_pair(index(GeometryTypes::hexahedron),  [&]() { return LagrangeCubeLocalFiniteElement<D,R,dim,compileTimeOrder>(runTimeOrder_); }),
            std::make_pair(prismIndex,                        [=]() { return LagrangePrismLocalFiniteElement<D,R,compileTimeOrder>(prismOrder); }),
            std::make_pair(pyramidIndex,                      [=]() { return LagrangePyramidLocalFiniteElement<D,R,compileTimeOrder>(pyramidOrder); })
          );
        }
        else
          return std::make_tuple(
            std::make_pair(index(GeometryTypes::simplex(dim)), [&]() { return LagrangeSimplexLocalFiniteElement<D,R,dim,compileTimeOrder>(runTimeOrder_); }),
            std::make_pair(index(GeometryTypes::cube(dim)),    [&]() { return LagrangeCubeLocalFiniteElement<D,R,dim,compileTimeOrder>(runTimeOrder_); })
          );
      }
    }
  };

} // namespace Impl



/** \brief A cache that stores all available Pk/Qk like local finite elements for the given dimension and order
 *
 * An interface for dealing with different vertex orders is currently missing.
 *
 * \tparam D Type used for domain coordinates
 * \tparam R Type used for shape function values
 * \tparam dim Element dimension
 * \tparam compileTimeOrder Polynomial order of the shape functions (or -1 for run-time order)
 *
 * The cached finite element implementations can be obtained using get(GeometryType).
 * If the \ref compileTimeOrder provided as template parameter is `-1`, the order
 * has to be provided as constructor argument.
 */
template<class D, class R, std::size_t dim, int compileTimeOrder = -1>
using LagrangeLocalFiniteElementCache = LocalFiniteElementVariantCache<Impl::ImplementedLagrangeFiniteElements<D,R,dim, compileTimeOrder>>;



} // namespace Dune




#endif // DUNE_LOCALFUNCTIONS_LAGRANGE_LAGRANGELFECACHE_HH
