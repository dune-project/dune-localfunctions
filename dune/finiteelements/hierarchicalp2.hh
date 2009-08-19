// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_HIERARCHICAL_P2_LOCALFINITEELEMENT_HH
#define DUNE_HIERARCHICAL_P2_LOCALFINITEELEMENT_HH

#include <dune/common/geometrytype.hh>

#include "common/localfiniteelement.hh"
#include "hierarchicalp2/hierarchicalp2/p2localbasis.hh"
#include "pk2d/pk2dlocalcoefficients.hh"
#include "pk2d/pk2dlocalinterpolation.hh"
#include "pk2d/pk2dlocalbasis.hh"
#include "pk3d/pk3dlocalcoefficients.hh"
#include "pk3d/pk3dlocalinterpolation.hh"
#include "pk3d/pk3dlocalbasis.hh"

namespace Dune
{

  /** \brief The Hierarchical P2 basis
   */
  template<class D, class R, int dim>
  class HierarchicalP2LocalFiniteElement
  {
  public:
    HierarchicalP2LocalFiniteElement()
    {
      DUNE_THROW(Dune::NotImplemented,"RefinedP1LocalFiniteElement not implemented for dim > 3.");
    }
  };

  /** \todo Please doc me !
   */
  template<class D, class R>
  class HierarchicalP2LocalFiniteElement<D,R,2> : LocalFiniteElementInterface<
                                                      LocalFiniteElementTraits<HierarchicalP2LocalBasis<D,R,2>,
                                                          Pk2DLocalCoefficients<2>,
                                                          Pk2DLocalInterpolation<Pk2DLocalBasis<D,R,2> > >
#ifndef DUNE_VIRTUAL_SHAPEFUNCTIONS
                                                      , HierarchicalP2LocalFiniteElement<D,R,2>
#endif
                                                      >
  {
  public:
    /** \todo Please doc me !
     */
    typedef LocalFiniteElementTraits<HierarchicalP2LocalBasis<D,R,2>,
        Pk2DLocalCoefficients<2>,
        Pk2DLocalInterpolation<Pk2DLocalBasis<D,R,2> > > Traits;

    /** \todo Please doc me !
     */
    HierarchicalP2LocalFiniteElement ()
    {
      gt.makeTriangle();
    }

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
    HierarchicalP2LocalBasis<D,R,2> basis;
    Pk2DLocalCoefficients<2> coefficients;
    Pk2DLocalInterpolation<Pk2DLocalBasis<D,R,2> > interpolation;
    GeometryType gt;
  };

  /** \todo Please doc me !
   */
  template<class D, class R>
  class HierarchicalP2LocalFiniteElement<D,R,3> : LocalFiniteElementInterface<
                                                      LocalFiniteElementTraits<HierarchicalP2LocalBasis<D,R,3>,
                                                          Pk3DLocalCoefficients<2>,
                                                          Pk3DLocalInterpolation<Pk3DLocalBasis<D,R,2> > >
#ifndef DUNE_VIRTUAL_SHAPEFUNCTIONS
                                                      , HierarchicalP2LocalFiniteElement<D,R,3>
#endif
                                                      >
  {
  public:
    /** \todo Please doc me !
     */
    typedef LocalFiniteElementTraits<HierarchicalP2LocalBasis<D,R,3>,
        Pk3DLocalCoefficients<2>,
        Pk3DLocalInterpolation<Pk3DLocalBasis<D,R,2> > > Traits;

    /** \todo Please doc me !
     */
    HierarchicalP2LocalFiniteElement ()
    {
      gt.makeTetrahedron();
    }

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
    HierarchicalP2LocalBasis<D,R,3> basis;
    Pk3DLocalCoefficients<2> coefficients;
    Pk3DLocalInterpolation<Pk3DLocalBasis<D,R,2> > interpolation;
    GeometryType gt;
  };

}

#endif
