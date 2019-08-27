// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_Q1_LOCALFINITEELEMENT_HH
#define DUNE_Q1_LOCALFINITEELEMENT_HH

#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include <dune/localfunctions/common/localtoglobaladaptors.hh>
#include <dune/localfunctions/lagrange/lagrangecube.hh>

namespace Dune
{

  /** \brief The local Q1 finite element on cubes
      \tparam D Domain data type
      \tparam R Range data type
      \tparam dim Dimension of the cube

      \deprecated This class is deprecated!  Please use LagrangeCubeLocalFiniteElement instead.
   */
  // The test test-q1.cc triggers compiler bugs in gcc-6 and earlier when Q1LocalFiniteElement
  // is simply redefined as LagrangeCubeLocalFiniteElement: it grabs more and more main memory,
  // and eventually stalls the machine.  Since I didn't find the real cause for this let's just
  // keep the relevant parts of Q1LocalFiniteElement until we can retire the problematic gcc versions.
#if !defined __GNUC__ || __GNUC__ > 6
  template<class D, class R, int dim>
  using Q1LocalFiniteElement = LagrangeCubeLocalFiniteElement<D,R,dim,1>;
#else
  template <int dim>
  class Q1LocalCoefficients
  {
  public:
    //! \brief Standard constructor
    Q1LocalCoefficients () : li(1<<dim)
    {
      for (std::size_t i=0; i<(1<<dim); i++)
        li[i] = LocalKey(i,dim,0);
    }

    //! number of coefficients
    std::size_t size () const
    {
      return 1<<dim;
    }

    //! get i'th index
    const LocalKey& localKey (std::size_t i) const
    {
      return li[i];
    }

  private:
    std::vector<LocalKey> li;
  };

  template<class D, class R, int dim>
  class Q1LocalFiniteElement
  {
  public:
    // user-defined default constructor is required for clang 3.8,
    // see https://gitlab.dune-project.org/core/dune-localfunctions/merge_requests/60
    /** default constructor */
    Q1LocalFiniteElement() {}

    /** \todo Please doc me !
     */
    typedef LocalFiniteElementTraits<Impl::LagrangeCubeLocalBasis<D,R,dim,1>,Q1LocalCoefficients<dim>,
        Impl::LagrangeCubeLocalInterpolation<Impl::LagrangeCubeLocalBasis<D,R,dim,1> > > Traits;

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
    static constexpr GeometryType type ()
    {
      return GeometryTypes::cube(dim);
    }

  private:
    Impl::LagrangeCubeLocalBasis<D,R,dim,1> basis;
    Q1LocalCoefficients<dim> coefficients;
    Impl::LagrangeCubeLocalInterpolation<Impl::LagrangeCubeLocalBasis<D,R,dim,1> > interpolation;
  };
#endif

  //! Factory for global-valued Q1 elements
  /**
   * \tparam Geometry Type of the geometry.  Used to extract the domain field
   *                  type and the dimension.
   * \tparam RF       Range field type.
   */
  template<class Geometry, class RF>
  class Q1FiniteElementFactory :
    public ScalarLocalToGlobalFiniteElementAdaptorFactory<
#if !defined __GNUC__ || __GNUC__ > 6
        LagrangeCubeLocalFiniteElement<
            typename Geometry::ctype, RF, Geometry::mydimension, 1
            >,
#else
        Q1LocalFiniteElement<
            typename Geometry::ctype, RF, Geometry::mydimension
            >,
#endif
        Geometry
        >
  {
#if !defined __GNUC__ || __GNUC__ > 6
    typedef LagrangeCubeLocalFiniteElement<
        typename Geometry::ctype, RF, Geometry::mydimension, 1
        > LFE;
#else
    typedef Q1LocalFiniteElement<
        typename Geometry::ctype, RF, Geometry::mydimension
        > LFE;
#endif
    typedef ScalarLocalToGlobalFiniteElementAdaptorFactory<LFE, Geometry> Base;

    static const LFE lfe;

  public:
    //! default constructor
    Q1FiniteElementFactory() : Base(lfe) {}
  };

  template<class Geometry, class RF>
  const typename Q1FiniteElementFactory<Geometry, RF>::LFE
  Q1FiniteElementFactory<Geometry, RF>::lfe;
}

#endif
