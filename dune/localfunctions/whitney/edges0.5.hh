// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

#ifndef DUNE_LOCALFUNCTIONS_WHITNEY_EDGES0_5_HH
#define DUNE_LOCALFUNCTIONS_WHITNEY_EDGES0_5_HH

#include <dune/geometry/type.hh>

#include <dune/localfunctions/whitney/edges0.5/basis.hh>
#include <dune/localfunctions/whitney/edges0.5/coefficients.hh>
#include <dune/localfunctions/whitney/edges0.5/interpolation.hh>

namespace Dune {

  //////////////////////////////////////////////////////////////////////
  //
  //  FiniteElement
  //

  //! FiniteElement for lowest order edge elements on simplices
  /**
   * Uses the representation
   * \f[
   *    \mathbf N^i=(L^{i_0}\nabla L^{i_1}-
   *                 L^{i_1}\nabla L^{i_0})\ell^i
   * \f]
   * where \f$L^k\f$ is the P1 shape function for vertex \f$k\f$, \f$i_0\f$
   * and \f$i_1\f$ are the indices of the vertices of edge \f$i\f$ and
   * \f$\ell^i\f$ is the length of edge \f$i\f$.
   *
   * \ingroup Whitney
   *
   * \tparam D   Type to represent the field in the domain.
   * \tparam R   Type to represent the field in the range.
   * \tparam dim Dimension of both domain and range.
   *
   * \nosubgrouping
   */
  template<class Geometry, class RF>
  class EdgeS0_5FiniteElement {
  public:
    /**
     * \implements FiniteElementInterface::Traits
     */
    struct Traits {
      typedef EdgeS0_5Basis<Geometry, RF> Basis;
      typedef EdgeS0_5Interpolation<Geometry,
          typename Basis::Traits> Interpolation;
      typedef EdgeS0_5Coefficients<Geometry::mydimension> Coefficients;
    };

  private:
    typename Traits::Basis basis_;
    typename Traits::Interpolation interpolation_;
    static const typename Traits::Coefficients& coefficients_;

  public:
    //! Constructor
    /**
     * \copydetails EdgeS0_5Basis::EdgeS0_5Basis(const Geometry& geo, const VertexOrder& vertexOrder)
     */
    template<class VertexOrder>
    EdgeS0_5FiniteElement(const Geometry& geo,
                          const VertexOrder& vertexOrder) :
      basis_(geo, vertexOrder), interpolation_(geo, vertexOrder)
    { }

    //! return reference to the basis object
    const typename Traits::Basis& basis() const { return basis_; }
    //! return reference to the interpolation object
    const typename Traits::Interpolation& interpolation() const
    { return interpolation_; }
    //! return reference to the coefficients object
    const typename Traits::Coefficients& coefficients() const
    { return coefficients_; }
    //! return geometry type of this element
    static constexpr GeometryType type() { return GeometryTypes::simplex(Geometry::mydimension); }
  };

  template<class Geometry, class RF>
  const typename EdgeS0_5FiniteElement<Geometry, RF>::Traits::Coefficients&
  EdgeS0_5FiniteElement<Geometry, RF>::coefficients_ =
    typename Traits::Coefficients();

  ////////////////////////////////////////////////////////////////////////
  //
  // Factory
  //

  //! Factory for EdgeS0_5FiniteElement objects
  /**
   * Constructs EdgeS0_5FiniteElement objects given a geometry and a vertex
   * ordering.
   *
   * \tparam Geometry Geometry for the local to global transformation.
   * \tparam RF       Field type of the range.
   *
   * \implements FiniteElementFactoryInterface
   */
  template<class Geometry, class RF>
  struct EdgeS0_5FiniteElementFactory {
    typedef EdgeS0_5FiniteElement<Geometry, RF> FiniteElement;

    //! construct the factory
    /**
     * \param geometry    The geometry object to use for adaption.
     * \param vertexOrder The global ordering of the vertices within the grid,
     *                    used to determine orientation of the edges.  This
     *                    vertexOrder object must support codim=0.
     *
     * \note The returned object stores the reference to the geometry passed
     *       here.  Any use of the returned value after this references has
     *       become invalid results in undefined behaviour.  The exception is
     *       that the destructor of this class may still be called.  The
     *       information contained in the vertexOrder object is extracted and
     *       the object is no longer needed after the contructor returns.  No
     *       reference to internal data of the factory is stored.
     */
    template<class VertexOrder>
    const FiniteElement make(const Geometry& geometry,
                             const VertexOrder& vertexOrder)
    { return FiniteElement(geometry, vertexOrder); }
  };

} // namespace Dune

#endif // DUNE_LOCALFUNCTIONS_WHITNEY_EDGES0_5_HH
