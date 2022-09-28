// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_PK1DLOCALFINITEELEMENT_HH
#define DUNE_PK1DLOCALFINITEELEMENT_HH

#include <cstddef>

#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include <dune/localfunctions/common/localtoglobaladaptors.hh>
#include <dune/localfunctions/lagrange/lagrangesimplex.hh>

#warning This header is deprecated

namespace Dune
{

  /** \brief Lagrange finite element on the unit interval with arbitrary compile-time order

      \deprecated This class is obsolete. Please use LagrangeSimplexLocalFiniteElement instead!
   */
  template<class D, class R, unsigned int k>
  using Pk1DLocalFiniteElement
    [[deprecated("use LagrangeSimplexLocalFiniteElement instead")]]
    = LagrangeSimplexLocalFiniteElement<D,R,1,k>;


  //! Langrange finite element of arbitrary order on triangles
  /**
   * \tparam Geometry Geometry for the local to global transformation.
   * \tparam RF       Field type of the range.
   * \tparam k        Maximum polynomial order of the base functions.
   *
   * \implements FiniteElementInterface
   */
  template<class Geometry, class RF, std::size_t k>
  class Pk1DFiniteElement {
    typedef typename Geometry::ctype DF;
    typedef Impl::LagrangeSimplexLocalBasis<DF,RF,1,k> LocalBasis;
    typedef Impl::LagrangeSimplexLocalInterpolation<LocalBasis> LocalInterpolation;

  public:
    /**
     * \implements FiniteElementInterface::Traits
     */
    struct Traits {
      typedef ScalarLocalToGlobalBasisAdaptor<LocalBasis, Geometry> Basis;
      typedef LocalToGlobalInterpolationAdaptor<
          LocalInterpolation,
          typename Basis::Traits
          > Interpolation;
      typedef Impl::LagrangeSimplexLocalCoefficients<1,k> Coefficients;
    };

  private:
    static const GeometryType gt;
    static const LocalBasis localBasis;
    static const LocalInterpolation localInterpolation;

    typename Traits::Basis basis_;
    typename Traits::Interpolation interpolation_;
    typename Traits::Coefficients coefficients_;

  public:
    //! construct a Pk1DFiniteElement
    /**
     * \param geometry    The geometry object to use for adaption.
     * \param vertexOrder The global ordering of the vertices within the grid,
     *                    used to determine orientation of the edges.  This
     *                    vertexOrder object must support codim=0.
     *
     * \note This class stores the reference to the geometry passed here.  Any
     *       use of this class after this references has become invalid
     *       results in undefined behaviour.  The exception is that the
     *       destructor of this class may still be called.  The information
     *       contained in the vertexOrder object is extracted and the object
     *       is no longer needed after the contructor returns.
     */
    template<class VertexOrder>
    Pk1DFiniteElement(const Geometry &geometry,
                      const VertexOrder& vertexOrder) :
      basis_(localBasis, geometry), interpolation_(localInterpolation),
      coefficients_(vertexOrder.begin(0, 0))
    { }

    const typename Traits::Basis& basis() const { return basis_; }
    const typename Traits::Interpolation& interpolation() const
    { return interpolation_; }
    const typename Traits::Coefficients& coefficients() const
    { return coefficients_; }
    const GeometryType &type() const { return gt; }
  };

  template<class Geometry, class RF, std::size_t k>
  const GeometryType
  Pk1DFiniteElement<Geometry, RF, k>::gt(GeometryTypes::simplex(2));

  template<class Geometry, class RF, std::size_t k>
  const typename Pk1DFiniteElement<Geometry, RF, k>::LocalBasis
  Pk1DFiniteElement<Geometry, RF, k>::localBasis = LocalBasis();

  template<class Geometry, class RF, std::size_t k>
  const typename Pk1DFiniteElement<Geometry, RF, k>::LocalInterpolation
  Pk1DFiniteElement<Geometry, RF, k>::localInterpolation =
    LocalInterpolation();

  //! Factory for Pk1DFiniteElement objects
  /**
   * Constructs Pk1DFiniteElement objects given a geometry and a vertex
   * ordering.
   *
   * \tparam Geometry Geometry for the local to global transformation.
   * \tparam RF       Field type of the range.
   * \tparam k        Maximum polynomial order of the base functions.
   *
   * \implements FiniteElementFactoryInterface
   */
  template<class Geometry, class RF, std::size_t k>
  struct Pk1DFiniteElementFactory {
    typedef Pk1DFiniteElement<Geometry, RF, k> FiniteElement;

    //! construct Pk1DFiniteElementFactory
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
     *       the object is no longer needed after the constructor returns.  No
     *       reference to internal data of the factory is stored.
     */
    template<class VertexOrder>
    const FiniteElement make(const Geometry& geometry,
                             const VertexOrder& vertexOrder)
    { return FiniteElement(geometry, vertexOrder); }
  };
}

#endif
