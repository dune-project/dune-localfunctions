// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

#ifndef DUNE_LOCALFUNCTIONS_MONOMIAL_HH
#define DUNE_LOCALFUNCTIONS_MONOMIAL_HH

#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <memory>
#include <vector>

#include <dune/geometry/type.hh>

#include "common/localfiniteelementtraits.hh"
#include "common/localtoglobaladaptors.hh"
#include "monomial/monomiallocalbasis.hh"
#include "monomial/monomiallocalcoefficients.hh"
#include "monomial/monomiallocalinterpolation.hh"

namespace Dune
{


  /** \brief Monomial basis for discontinuous Galerkin methods
   *
   * Be careful: Although MonomialLocalInterpolation::interpolate
   * uses an L^2 projection it is unstable for higher polynomial degrees.
   *
   * \ingroup Monomial
   *
   * \tparam D Type used for coordinates
   * \tparam R Type used for shape function values
   * \tparam d Dimension of the element
   * \tparam p Order of the basis
   */
  template<class D, class R, int d, int p>
  class MonomialLocalFiniteElement
  {
    constexpr static int static_size = MonomialLocalBasis<D,R,d,p>::size();

  public:
    /** Traits class
     */
    typedef LocalFiniteElementTraits<
        MonomialLocalBasis<D,R,d,p>,
        MonomialLocalCoefficients<static_size>,
        MonomialLocalInterpolation<MonomialLocalBasis<D,R,d,p>,static_size>
        > Traits;

    //! Construct a MonomLocalFiniteElement
    MonomialLocalFiniteElement (const GeometryType &gt_)
      : basis(), interpolation(gt_, basis), gt(gt_)
    {}

    /** \todo Please doc me !
     */
    const typename Traits::LocalBasisType& localBasis () const
    {
      return basis;
    }

    /** \todo Please doc me !
     */
    const typename Traits::LocalCoefficientsType& localCoefficients () const
    {
      return coefficients;
    }

    /** \todo Please doc me !
     */
    const typename Traits::LocalInterpolationType& localInterpolation () const
    {
      return interpolation;
    }

    /** \brief Number of shape functions in this finite element */
    unsigned int size () const
    {
      return basis.size();
    }

    /** \todo Please doc me !
     */
    GeometryType type () const
    {
      return gt;
    }

  private:
    MonomialLocalBasis<D,R,d,p> basis;
    MonomialLocalCoefficients<static_size> coefficients;
    MonomialLocalInterpolation<MonomialLocalBasis<D,R,d,p>,static_size> interpolation;
    GeometryType gt;
  };

  //! Factory for global-valued MonomFiniteElement objects
  /**
   * Constructs MonomialFiniteElement objects given a geometry.
   *
   * \tparam Geometry Geometry for the local to global transformation.
   * \tparam RF       Field type of the range.
   * \tparam p        Order of the basis.
   *
   * \implements FiniteElementFactoryInterface
   *
   * \note There is no real MonomFiniteElement, only the FiniteElement typedef
   *       inside this class.
   */
  template<class Geometry, class RF, std::size_t p>
  class MonomialFiniteElementFactory {
    typedef typename Geometry::ctype DF;
    static const std::size_t dim = Geometry::mydimension;

    typedef MonomialLocalFiniteElement<DF, RF, dim, p> LocalFE;

    std::vector<std::shared_ptr<const LocalFE> > localFEs;

    void init(const GeometryType &gt) {
      std::size_t index = gt.id() >> 1;
      if(localFEs.size() <= index)
        localFEs.resize(index+1);
      localFEs[index].reset(new LocalFE(gt));
    }

  public:
    typedef ScalarLocalToGlobalFiniteElementAdaptor<LocalFE, Geometry>
    FiniteElement;

    //! construct a MonomialFiniteElementFactory from a list of GeometryType's
    /**
     * \param begin Begin of a range of geometry types.
     * \param end   End of a range of geometry types.
     */
    template<class ForwardIterator>
    MonomialFiniteElementFactory(const ForwardIterator &begin,
                              const ForwardIterator &end)
    {
      for(ForwardIterator it = begin; it != end; ++it)
        init(*it);
    }

    //! construct a MonomialFiniteElementFactory from a single GeometryType
    /**
     * \param gt GeometryType to construct elements with
     */
    MonomialFiniteElementFactory(const GeometryType &gt)
    { init(gt); }

    //! construct a MonomFiniteElementFactory for all applicable GeometryType's
    /**
     * \note This constructor only works for dimensions up to and including 3.
     */
    MonomialFiniteElementFactory() {
      static_assert(dim <= 3, "MonomFiniteElementFactory knows the "
                    "available geometry types only up to dimension 3");

      GeometryType gt;
      switch(dim) {
      case 0 :
        gt = Dune::GeometryTypes::vertex;        init(gt);
        break;
      case 1 :
        gt = Dune::GeometryTypes::line;          init(gt);
        break;
      case 2 :
        gt = Dune::GeometryTypes::triangle;      init(gt);
        gt = Dune::GeometryTypes::quadrilateral; init(gt);
        break;
      case 3 :
        gt = Dune::GeometryTypes::tetrahedron;   init(gt);
        gt = Dune::GeometryTypes::pyramid;       init(gt);
        gt = Dune::GeometryTypes::prism;         init(gt);
        gt = Dune::GeometryTypes::hexahedron;    init(gt);
        break;
      default :
        // this should never happen -- it should be caught by the static
        // assert above.
        std::abort();
      };
    }

    //! construct a global-valued MonomFiniteElement
    /**
     * \param geometry The geometry object to use for adaption.
     *
     * \note The returned object stores the reference to the geometry passed
     *       here as well as references to internal data of this factory.  Any
     *       use of the returned value after the geometry reference or the
     *       factory object was become invalid results in undefined behaviour.
     *       The exception is that the destructor of the returned value may
     *       still be called.
     */
    const FiniteElement make(const Geometry& geometry) {
      std::size_t index = geometry.type().id() >> 1;
      assert(localFEs.size() > index && localFEs[index]);
      return FiniteElement(*localFEs[index], geometry);
    }
  };
}

#endif // DUNE_LOCALFUNCTIONS_MONOMIAL_HH
