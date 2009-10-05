// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_Q1_LOCALFINITEELEMENT_HH
#define DUNE_Q1_LOCALFINITEELEMENT_HH

#include <dune/finiteelements/common/localfiniteelement.hh>
#include <dune/finiteelements/q1/q1localbasis.hh>
#include <dune/finiteelements/q1/q1localcoefficients.hh>
#include <dune/finiteelements/q1/q1localinterpolation.hh>

namespace Dune
{

  /** \todo Please doc me !
   */
  template<class D, class R, int dim>
  class Q1LocalFiniteElement
    : public LocalFiniteElementInterface<
          LocalFiniteElementTraits<Q1LocalBasis<D,R,dim>,Q1LocalCoefficients<dim>,
              Q1LocalInterpolation<dim,Q1LocalBasis<D,R,dim> > >
#ifndef DUNE_VIRTUAL_SHAPEFUNCTIONS
          ,Q12DLocalFiniteElement<D,R,dim>
#endif
          >
  {
  public:
    /** \todo Please doc me !
     */
    typedef LocalFiniteElementTraits<Q1LocalBasis<D,R,dim>,Q1LocalCoefficients<dim>,
        Q1LocalInterpolation<dim,Q1LocalBasis<D,R,dim> > > Traits;

    /** \todo Please doc me !
     */
    Q1LocalFiniteElement ()
    {
      gt.makeCube(dim);
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
    Q1LocalBasis<D,R,dim> basis;
    Q1LocalCoefficients<dim> coefficients;
    Q1LocalInterpolation<dim,Q1LocalBasis<D,R,dim> > interpolation;
    GeometryType gt;
  };

}

#endif
