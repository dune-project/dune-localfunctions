// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_Q2_LOCALFINITEELEMENT_HH
#define DUNE_Q2_LOCALFINITEELEMENT_HH

#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include <dune/localfunctions/common/localtoglobaladaptors.hh>
#include "q2/q2localbasis.hh"
#include "q2/q2localcoefficients.hh"
#include "q2/q2localinterpolation.hh"

namespace Dune
{

  /** \brief 2nd-order Lagrangian finite elements on hypercubes
   * \tparam D Type used for coordinates
   * \tparam R Type used for function values
   * \tparam dim Dimension of the reference cube
   */
  template<class D, class R, int dim>
  class Q2LocalFiniteElement
  {
  public:
    /** \todo Please doc me !
     */
    typedef LocalFiniteElementTraits<Q2LocalBasis<D,R,dim>,Q2LocalCoefficients<dim>,
        Q2LocalInterpolation<Q2LocalBasis<D,R,dim> > > Traits;

    /** \brief Default constructor
     */
    Q2LocalFiniteElement ()
    {
      gt.makeCube(dim);
    }

    /** \brief Get the actual shape functions
     */
    const typename Traits::LocalBasisType& localBasis () const
    {
      return basis;
    }

    /** \todo Please doc me!
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

    /** \brief Get the type of the reference element
     */
    GeometryType type () const
    {
      return gt;
    }

  private:
    Q2LocalBasis<D,R,dim> basis;
    Q2LocalCoefficients<dim> coefficients;
    Q2LocalInterpolation<Q2LocalBasis<D,R,dim> > interpolation;
    GeometryType gt;
  };

  //! Factory for global-valued Q23D elements
  /**
   * \tparam Geometry Type of the geometry.  Used to extract the domain field
   *                  type.
   * \tparam RF       Range field type.
   */
  template<class Geometry, class RF>
  class Q2FiniteElementFactory :
    public ScalarLocalToGlobalFiniteElementAdaptorFactory<
        Q2LocalFiniteElement<typename Geometry::ctype, RF, Geometry::mydimension>, Geometry
        >
  {
    typedef Q2LocalFiniteElement<typename Geometry::ctype, RF, Geometry::mydimension> LFE;
    typedef ScalarLocalToGlobalFiniteElementAdaptorFactory<LFE, Geometry> Base;

    static const LFE lfe;

  public:
    //! default constructor
    Q2FiniteElementFactory() : Base(lfe) {}
  };

  template<class Geometry, class RF>
  const typename Q2FiniteElementFactory<Geometry, RF>::LFE
  Q2FiniteElementFactory<Geometry, RF>::lfe;
}

#endif
