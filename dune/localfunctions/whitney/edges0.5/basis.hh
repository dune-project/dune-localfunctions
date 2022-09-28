// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

#ifndef DUNE_LOCALFUNCTIONS_WHITNEY_EDGES0_5_BASIS_HH
#define DUNE_LOCALFUNCTIONS_WHITNEY_EDGES0_5_BASIS_HH

#include <cstddef>
#include <vector>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <dune/localfunctions/common/localtoglobaladaptors.hh>
#include <dune/localfunctions/lagrange/lagrangesimplex.hh>
#include <dune/localfunctions/whitney/edges0.5/common.hh>

namespace Dune {

  //////////////////////////////////////////////////////////////////////
  //
  //  Basis
  //

  //! Basis for order 0.5 (lowest order) edge elements on simplices
  /**
   * @ingroup BasisImplementation
   *
   * \tparam Geometry Type of the local-to-global map.
   * \tparam RF       Type to represent the field in the range.
   *
   * \nosubgrouping
   */
  template<class Geometry, class RF>
  class EdgeS0_5Basis :
    private EdgeS0_5Common<Geometry::mydimension, typename Geometry::ctype>
  {
  public:
    //! \brief export type traits for function signature
    struct Traits {
      typedef typename Geometry::ctype DomainField;
      static const std::size_t dimDomainLocal = Geometry::mydimension;
      static const std::size_t dimDomainGlobal = Geometry::coorddimension;
      typedef FieldVector<DomainField, dimDomainLocal> DomainLocal;
      typedef FieldVector<DomainField, dimDomainGlobal> DomainGlobal;

      typedef RF RangeField;
      static const std::size_t dimRange = dimDomainLocal;
      typedef FieldVector<RangeField, dimRange> Range;

      typedef FieldMatrix<RangeField, dimRange, dimDomainGlobal> Jacobian;
    };

  private:
    typedef Dune::Impl::LagrangeSimplexLocalBasis<typename Traits::DomainField,
        typename Traits::RangeField,
        Traits::dimDomainLocal,
        1    // Polynomial order
        > P1LocalBasis;
    typedef ScalarLocalToGlobalBasisAdaptor<P1LocalBasis, Geometry> P1Basis;

    static const P1LocalBasis& p1LocalBasis;
    static const std::size_t dim = Traits::dimDomainLocal;

    typedef EdgeS0_5Common<dim, typename Geometry::ctype> Base;
    using Base::refelem;
    using Base::s;

    // global values of the Jacobians (gradients) of the p1 basis
    std::vector<typename P1Basis::Traits::Jacobian> p1j;
    // edge sizes and orientations
    std::vector<typename Traits::DomainField> edgel;

  public:
    //! Construct an EdgeS0_5Basis
    /**
     * \param geo         Geometry of the element to contruct a local basis
     *                    for.
     * \param vertexOrder Vertex ordering information.  Only the vertex order
     *                    on the dim=1 sub-entities (edges) is required.
     */
    template<typename VertexOrder>
    EdgeS0_5Basis(const Geometry& geo, const VertexOrder& vertexOrder) :
      p1j(s, typename P1Basis::Traits::Jacobian(0)), edgel(s)
    {
      // use some arbitrary position to evaluate jacobians, they are constant
      static const typename Traits::DomainLocal xl(0);

      // precompute Jacobian (gradients) of the p1 element
      P1Basis(p1LocalBasis, geo).evaluateJacobian(xl, p1j);

      // calculate edge sizes and orientations
      for(std::size_t i = 0; i < s; ++i) {
        edgel[i] = (geo.corner(refelem.subEntity(i,dim-1,0,dim))-
                    geo.corner(refelem.subEntity(i,dim-1,1,dim))
                    ).two_norm();
        const typename VertexOrder::iterator& edgeVertexOrder =
          vertexOrder.begin(dim-1, i);
        if(edgeVertexOrder[0] > edgeVertexOrder[1])
          edgel[i] *= -1;
      }
    }

    //! number of shape functions
    std::size_t size () const { return s; }

    //! Evaluate all shape functions
    void evaluateFunction(const typename Traits::DomainLocal& xl,
                          std::vector<typename Traits::Range>& out) const
    {
      out.assign(s, typename Traits::Range(0));

      // compute p1 values -- use the local basis directly for that, local and
      // global values are identical for scalars
      std::vector<typename P1LocalBasis::Traits::RangeType> p1v;
      p1LocalBasis.evaluateFunction(xl, p1v);

      for(std::size_t i = 0; i < s; i++) {
        const std::size_t i0 = refelem.subEntity(i,dim-1,0,dim);
        const std::size_t i1 = refelem.subEntity(i,dim-1,1,dim);
        out[i].axpy( p1v[i0], p1j[i1][0]);
        out[i].axpy(-p1v[i1], p1j[i0][0]);
        out[i] *= edgel[i];
      }
    }

    //! Evaluate all Jacobians
    void evaluateJacobian(const typename Traits::DomainLocal&,
                          std::vector<typename Traits::Jacobian>& out) const
    {
      out.resize(s);

      for(std::size_t i = 0; i < s; i++) {
        const std::size_t i0 = refelem.subEntity(i,dim-1,0,dim);
        const std::size_t i1 = refelem.subEntity(i,dim-1,1,dim);
        for(std::size_t j = 0; j < dim; j++)
          for(std::size_t k = 0; k < dim; k++)
            out[i][j][k] = edgel[i] *
                           (p1j[i0][0][k]*p1j[i1][0][j]-p1j[i1][0][k]*p1j[i0][0][j]);
      }
    }

    //! \brief Evaluate partial derivatives of all shape functions
    void partial (const std::array<unsigned int, dim>& order,
                  const typename Traits::DomainLocal& in,         // position
                  std::vector<typename Traits::Range>& out) const      // return value
    {
      auto totalOrder = std::accumulate(order.begin(), order.end(), 0);
      if (totalOrder == 0) {
        evaluateFunction(in, out);
      } else if (totalOrder==1) {
        auto const k = std::distance(order.begin(), std::find(order.begin(), order.end(), 1));
        out.resize(size());

        for (std::size_t i = 0; i < s; i++)
        {
          const std::size_t i0 = refelem.subEntity(i,dim-1,0,dim);
          const std::size_t i1 = refelem.subEntity(i,dim-1,1,dim);
          for(std::size_t j = 0; j < dim; j++)
            out[i][j] = edgel[i] *
              (p1j[i0][0][k]*p1j[i1][0][j] - p1j[i1][0][k]*p1j[i0][0][j]);
        }
      } else {
        DUNE_THROW(NotImplemented, "Desired derivative order is not implemented");
      }
    }

    //! Polynomial order of the shape functions
    std::size_t order () const { return 1; }
  };

  template<class Geometry, class RF>
  const typename EdgeS0_5Basis<Geometry, RF>::P1LocalBasis&
  EdgeS0_5Basis<Geometry, RF>::p1LocalBasis = P1LocalBasis();

} // namespace Dune

#endif // DUNE_LOCALFUNCTIONS_WHITNEY_EDGES0_5_BASIS_HH
