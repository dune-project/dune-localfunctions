// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LAGRANGEFINITEELEMENT_HH
#define DUNE_LAGRANGEFINITEELEMENT_HH

#include <dune/common/geometrytype.hh>
#include <dune/grid/genericgeometry/conversion.hh>

#include <dune/finiteelements/common/localfiniteelement.hh>
#include <dune/finiteelements/lagrangebasis/lagrangepoints.hh>
#include <dune/finiteelements/lagrangebasis/lagrangebasis.hh>

namespace Dune
{

  template< unsigned int dimDomain, class D, class R, class CF=D >
  class LagrangeLocalFiniteElement
    : LocalFiniteElementInterface<
          LocalFiniteElementTraits< GenericLocalBasis<dimDomain,D,R,typename LagrangeBasisProvider< dimDomain, D, CF >::Basis > ,
              typename LagrangePointsCreator< D, dimDomain >::LocalCoefficients,
              typename LocalLagrangeInterpolationCreator< LagrangePointsCreator< D, dimDomain > >::LocalInterpolation >,
          LagrangeLocalFiniteElement<dimDomain,D,R,CF> >
  {
    typedef LagrangeBasisProvider< dimDomain, D, CF > BasisCreator;
    typedef Dune::LagrangePointsCreator< D, dimDomain > LagrangePointsCreator;
    typedef LagrangePointsCreator LocalCoefficientsCreator;
    typedef LocalLagrangeInterpolationCreator< LagrangePointsCreator > LocalInterpolationCreator;
    typedef FiniteElementProvider<BasisCreator,LocalCoefficientsCreator,LocalInterpolationCreator> FECreator;
    typedef typename FECreator::FiniteElement FiniteElement;
  public:
    /** \todo Please doc me !
     */
    typedef LocalFiniteElementTraits< GenericLocalBasis<dimDomain,D,R,typename LagrangeBasisProvider< dimDomain, D, CF >::Basis > ,
        typename LagrangePointsCreator::LocalCoefficients,
        typename LocalInterpolationCreator::LocalInterpolation > Traits;

    /** \todo Please doc me !
     */
    LagrangeLocalFiniteElement ( unsigned int topologyId,
                                 unsigned int order )
      : topologyId_(topologyId),
        order_(order),
        finiteElement_( FECreator::finiteElement(topologyId,order) ),
        localBasis_(finiteElement_.basis())
    {}

    /** \todo Please doc me !
     */
    const typename Traits::LocalBasisType& localBasis () const
    {
      return localBasis_;
    }

    /** \todo Please doc me !
     */
    const typename Traits::LocalCoefficientsType& localCoefficients () const
    {
      return finiteElement_.coefficients();
    }

    /** \todo Please doc me !
     */
    const typename Traits::LocalInterpolationType& localInterpolation () const
    {
      return finiteElement_.interpolation();
    }

    /** \todo Please doc me !
     */
    GeometryType type () const
    {
      if ( GenericGeometry::hasGeometryType( topologyId_, dimDomain ) )
        return GenericGeometry::geometryType( topologyId_, dimDomain );
      return GeometryType();
    }

    /** \todo Please doc me !
     */
    unsigned int topologyId () const
    {
      return topologyId_;
    }
  private:
    unsigned int topologyId_;
    unsigned int order_;
    const FiniteElement &finiteElement_;
    typename Traits::LocalBasisType localBasis_;
  };

}

#endif
