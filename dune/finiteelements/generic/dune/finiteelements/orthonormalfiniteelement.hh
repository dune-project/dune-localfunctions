// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ORTHONORMALFINITEELEMENT_HH
#define DUNE_ORTHONORMALFINITEELEMENT_HH

#include <dune/finiteelements/generic/localfiniteelement.hh>
#include <dune/finiteelements/generic/dglocalcoefficients.hh>
#include <dune/finiteelements/orthonormalbasis/l2interpolation.hh>
#include <dune/finiteelements/orthonormalbasis/orthonormalbasis.hh>

namespace Dune
{

  template< unsigned int dimDomain, class D, class R, class CF=D >
  class OrthonormalLocalFiniteElement
    : public GenericLocalFiniteElement< OrthonormalBasisProvider< dimDomain, D, CF >,
          DGLocalCoefficientsCreator< OrthonormalBasisProvider< dimDomain, D, CF > >,
          LocalL2InterpolationCreator< OrthonormalBasisProvider< dimDomain, D, CF > >,
          dimDomain,D,R,CF>
  {
    typedef GenericLocalFiniteElement< OrthonormalBasisProvider< dimDomain, R, CF >,
        DGLocalCoefficientsCreator< OrthonormalBasisProvider< dimDomain, R, CF > >,
        LocalL2InterpolationCreator< OrthonormalBasisProvider< dimDomain, R, CF > >,
        dimDomain,D,R,CF> Base;
    using Base::FECreator;
    using Base::FiniteElement;
  public:
    using Base::Traits;

    /** \todo Please doc me !
     */
    OrthonormalLocalFiniteElement ( unsigned int topologyId,
                                    unsigned int order )
      : Base(topologyId,order)
    {}
  };

}

#endif
