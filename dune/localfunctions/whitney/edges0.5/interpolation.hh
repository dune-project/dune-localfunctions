// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

#ifndef DUNE_LOCALFUNCTIONS_WHITNEY_EDGES0_5_INTERPOLATION_HH
#define DUNE_LOCALFUNCTIONS_WHITNEY_EDGES0_5_INTERPOLATION_HH

#include <cstddef>
#include <vector>

#include <dune/localfunctions/whitney/edges0.5/common.hh>
#include <dune/localfunctions/common/localinterpolation.hh>

namespace Dune {

  //////////////////////////////////////////////////////////////////////
  //
  // Interpolation
  //

  //! Interpolation for lowest order edge elements on simplices
  /**
   * \tparam Geometry Type of the local-to-global map.
   * \tparam RF       Type to represent the field in the range.
   *
   * \nosubgrouping
   */
  template<class Geometry, class Traits_>
  class EdgeS0_5Interpolation :
    private EdgeS0_5Common<Traits_::dimDomainLocal,
        typename Traits_::DomainField>
  {
  public:
    typedef Traits_ Traits;

  private:
    static const std::size_t dim = Traits::dimDomainLocal;
    typedef EdgeS0_5Common<dim, typename Traits::DomainField> Base;
    using Base::refelem;
    using Base::s;

    std::vector<typename Traits::DomainGlobal> edgev;

  public:
    //! constructor
    /**
     * \param geo         Geometry of the element to contruct a local basis
     *                    for.
     * \param vertexOrder Vertex ordering information.  Only the vertex order
     *                    on the dim=1 sub-entities (edges) is required.
     */
    template<typename VertexOrder>
    EdgeS0_5Interpolation(const Geometry& geo,
                          const VertexOrder& vertexOrder) :
      edgev(s)
    {
      for(std::size_t i = 0; i < s; ++i) {
        const std::size_t i0 = refelem.subEntity(i,dim-1,0,dim);
        const std::size_t i1 = refelem.subEntity(i,dim-1,1,dim);

        edgev[i] = geo.corner(i1);
        edgev[i] -= geo.corner(i0);
        edgev[i] /= edgev[i].two_norm();

        const typename VertexOrder::iterator& edgeVertexOrder =
          vertexOrder.begin(dim-1, i);
        if(edgeVertexOrder[0] > edgeVertexOrder[1])
          edgev[i] *= -1;
      }
    }

    //! Interpolation of a function
    template<typename F, typename C>
    void interpolate(const F& ff, std::vector<C>& out) const {
      typename Traits::Range y;

      auto&& f = Impl::makeFunctionWithCallOperator<std::decay_t<decltype(refelem.position(0,dim-1))>>(ff);

      out.resize(s);

      for(std::size_t i = 0; i < s; ++i) {
        y = f(refelem.position(i,dim-1));

        out[i] = y * edgev[i];
      }
    }
  };

} // namespace Dune

#endif // DUNE_LOCALFUNCTIONS_WHITNEY_EDGES0_5_INTERPOLATION_HH
