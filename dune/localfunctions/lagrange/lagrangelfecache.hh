// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LOCALFUNCTIONS_LAGRANGE_LAGRANGELFECACHE_HH
#define DUNE_LOCALFUNCTIONS_LAGRANGE_LAGRANGELFECACHE_HH

#include <tuple>
#include <utility>

#include <dune/geometry/type.hh>
#include <dune/geometry/typeindex.hh>

#include <dune/localfunctions/lagrange/p0.hh>
#include <dune/localfunctions/lagrange/pk.hh>
#include <dune/localfunctions/lagrange/qk.hh>
#include <dune/localfunctions/lagrange/prismp1.hh>
#include <dune/localfunctions/lagrange/prismp2.hh>
#include <dune/localfunctions/lagrange/pyramidp1.hh>
#include <dune/localfunctions/lagrange/pyramidp2.hh>
#include <dune/localfunctions/common/localfiniteelementvariantcache.hh>


namespace Dune {



namespace Impl {

  // Provide implemented Lagrange local finite elements

  template<class D, class R, std::size_t dim, std::size_t order>
  struct ImplementedLagrangeFiniteElements : public FixedDimLocalGeometryTypeIndex<dim>
  {
    using FixedDimLocalGeometryTypeIndex<dim>::index;
    static auto getImplementations()
    {
      return std::make_tuple(
        std::make_pair(index(GeometryTypes::simplex(dim)), []() { return PkLocalFiniteElement<D,R,dim,order>(); }),
        std::make_pair(index(GeometryTypes::cube(dim)),    []() { return QkLocalFiniteElement<D,R,dim,order>(); })
      );
    }
  };

  template<class D, class R, std::size_t dim>
  struct ImplementedLagrangeFiniteElements<D,R,dim,0> : public FixedDimLocalGeometryTypeIndex<dim>
  {
    using FixedDimLocalGeometryTypeIndex<dim>::index;
    static auto getImplementations()
    {
      return std::make_tuple(
        std::make_pair(index(GeometryTypes::simplex(dim)), []() { return P0LocalFiniteElement<D,R,dim>(GeometryTypes::simplex(dim)); }),
        std::make_pair(index(GeometryTypes::cube(dim)),    []() { return P0LocalFiniteElement<D,R,dim>(GeometryTypes::cube(dim)); }),
        std::make_pair(index(GeometryTypes::none(dim)),    []() { return P0LocalFiniteElement<D,R,dim>(GeometryTypes::none(dim)); })
      );
    }
  };

  template<class D, class R>
  struct ImplementedLagrangeFiniteElements<D,R,3,0> : public FixedDimLocalGeometryTypeIndex<3>
  {
    using FixedDimLocalGeometryTypeIndex<3>::index;
    static auto getImplementations()
    {
      return std::make_tuple(
        std::make_pair(index(GeometryTypes::tetrahedron), []() { return P0LocalFiniteElement<D,R,3>(GeometryTypes::tetrahedron); }),
        std::make_pair(index(GeometryTypes::hexahedron),  []() { return P0LocalFiniteElement<D,R,3>(GeometryTypes::hexahedron); }),
        std::make_pair(index(GeometryTypes::prism),       []() { return P0LocalFiniteElement<D,R,3>(GeometryTypes::prism); }),
        std::make_pair(index(GeometryTypes::pyramid),     []() { return P0LocalFiniteElement<D,R,3>(GeometryTypes::pyramid); })
      );
    }
  };

  template<class D, class R>
  struct ImplementedLagrangeFiniteElements<D,R,3,1> : public FixedDimLocalGeometryTypeIndex<3>
  {
    using FixedDimLocalGeometryTypeIndex<3>::index;
    static auto getImplementations()
    {
      return std::make_tuple(
        std::make_pair(index(GeometryTypes::tetrahedron), []() { return PkLocalFiniteElement<D,R,3,1>(); }),
        std::make_pair(index(GeometryTypes::hexahedron),  []() { return QkLocalFiniteElement<D,R,3,1>(); }),
        std::make_pair(index(GeometryTypes::prism),       []() { return PrismP1LocalFiniteElement<D,R>(); }),
        std::make_pair(index(GeometryTypes::pyramid),     []() { return PyramidP1LocalFiniteElement<D,R>(); })
      );
    }
  };

  template<class D, class R>
  struct ImplementedLagrangeFiniteElements<D,R,3,2> : public FixedDimLocalGeometryTypeIndex<3>
  {
    using FixedDimLocalGeometryTypeIndex<3>::index;
    static auto getImplementations()
    {
      return std::make_tuple(
        std::make_pair(index(GeometryTypes::tetrahedron), []() { return PkLocalFiniteElement<D,R,3,2>(); }),
        std::make_pair(index(GeometryTypes::hexahedron),  []() { return QkLocalFiniteElement<D,R,3,2>(); }),
        std::make_pair(index(GeometryTypes::prism),       []() { return PrismP2LocalFiniteElement<D,R>(); }),
        std::make_pair(index(GeometryTypes::pyramid),     []() { return PyramidP2LocalFiniteElement<D,R>(); })
      );
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
 * \tparam order Element order
 *
 * The cached finite element implementations can be obtained using get(GeometryType).
 */
template<class D, class R, std::size_t dim, std::size_t order>
using LagrangeLocalFiniteElementCache = LocalFiniteElementVariantCache<Impl::ImplementedLagrangeFiniteElements<D,R,dim,order>>;



} // namespace Dune




#endif // DUNE_LOCALFUNCTIONS_LAGRANGE_LAGRANGELFECACHE_HH
