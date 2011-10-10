// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_MONOMLOCALFINITEELEMENT_HH
#define DUNE_MONOMLOCALFINITEELEMENT_HH

#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <vector>

#include <dune/common/deprecated.hh>
#include <dune/common/shared_ptr.hh>
#include <dune/common/static_assert.hh>

#include <dune/geometry/type.hh>

#include "common/localfiniteelementtraits.hh"
#include "common/localtoglobaladaptors.hh"
#include "monom/monomlocalbasis.hh"
#include "monom/monomlocalcoefficients.hh"
#include "monom/monomlocalinterpolation.hh"

namespace Dune
{

  /** Monom basis for discontinuous Galerkin

     \tparam D Type used for coordinates
     \tparam R Type used for shape function values
     \tparam d Dimension of the element
     \tparam p Order of the basis
   */
  template<class D, class R, int d, int p>
  class MonomLocalFiniteElement
  {
    enum { static_size = MonomImp::Size<d,p>::val };

  public:
    /** Traits class
     */
    typedef LocalFiniteElementTraits<
        MonomLocalBasis<D,R,d,p>,
        MonomLocalCoefficients<static_size>,
        MonomLocalInterpolation<MonomLocalBasis<D,R,d,p>,static_size>
        > Traits;

    /** \todo Please doc me !
     */
    MonomLocalFiniteElement (GeometryType::BasicType basicType) DUNE_DEPRECATED
      : basis(), interpolation(basicType, basis), gt(basicType,d)
    {}

    //! Construct a MonomLocalFiniteElement
    MonomLocalFiniteElement (const GeometryType &gt_)
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

    /** \todo Please doc me !
     */
    GeometryType type () const
    {
      return gt;
    }

  private:
    MonomLocalBasis<D,R,d,p> basis;
    MonomLocalCoefficients<static_size> coefficients;
    MonomLocalInterpolation<MonomLocalBasis<D,R,d,p>,static_size> interpolation;
    GeometryType gt;
  };

  //! Factory for global-valued MonomFiniteElement objects
  /**
   * Constructs MonomFiniteElement objects given a geometry.
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
  class MonomFiniteElementFactory {
    typedef typename Geometry::ctype DF;
    static const std::size_t dim = Geometry::mydimension;

    typedef MonomLocalFiniteElement<DF, RF, dim, p> LocalFE;

    std::vector<shared_ptr<const LocalFE> > localFEs;

    void init(const GeometryType &gt) {
      std::size_t index = gt.id() >> 1;
      if(localFEs.size() <= index)
        localFEs.resize(index+1);
      localFEs[index].reset(new LocalFE(gt));
    }

  public:
    typedef ScalarLocalToGlobalFiniteElementAdaptor<LocalFE, Geometry>
    FiniteElement;

    //! construct a MonomFiniteElementFactory from a list of GeometryType's
    /**
     * \param begin Begin of a range of geometry types.
     * \param end   End of a range of geometry types.
     */
    template<class ForwardIterator>
    MonomFiniteElementFactory(const ForwardIterator &begin,
                              const ForwardIterator &end)
    {
      for(ForwardIterator it = begin; it != end; ++it)
        init(*it);
    }

    //! construct a MonomFiniteElementFactory from a single GeometryType
    /**
     * \param gt GeometryType to construct elements with
     */
    MonomFiniteElementFactory(const GeometryType &gt)
    { init(gt); }

    //! construct a MonomFiniteElementFactory for all applicable GeometryType's
    /**
     * \note This constructor only works for dimensions up to and including 3.
     */
    MonomFiniteElementFactory() {
      dune_static_assert(dim <= 3, "MonomFiniteElementFactory knows the "
                         "available geometry types only up to dimension 3");

      GeometryType gt;
      switch(dim) {
      case 0 :
        gt.makeVertex();        init(gt);
        break;
      case 1 :
        gt.makeLine();          init(gt);
        break;
      case 2 :
        gt.makeTriangle();      init(gt);
        gt.makeQuadrilateral(); init(gt);
        break;
      case 3 :
        gt.makeTetrahedron();   init(gt);
        gt.makePyramid();       init(gt);
        gt.makePrism();         init(gt);
        gt.makeHexahedron();    init(gt);
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

#endif // DUNE_MONOMLOCALFINITEELEMENT_HH
