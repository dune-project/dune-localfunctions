// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ORTHONORMALFINITEELEMENT_HH
#define DUNE_ORTHONORMALFINITEELEMENT_HH

#include <dune/finiteelements/generic/localfiniteelement.hh>
#include <dune/finiteelements/generic/dglocalcoefficients.hh>
#include <dune/finiteelements/generic/l2interpolation.hh>
#include <dune/finiteelements/orthonormalbasis/orthonormalbasis.hh>

namespace Dune
{

  template< unsigned int dimDomain, class D, class R,
      class SF=R, class CF=SF >

  class OrthonormalLocalFiniteElement
    : public GenericLocalFiniteElement< OrthonormalBasisFactory< dimDomain, SF, CF >,
          DGLocalCoefficientsFactory< OrthonormalBasisFactory< dimDomain, SF, CF > >,
          LocalL2InterpolationFactory< OrthonormalBasisFactory< dimDomain, SF, CF >,true >,
          dimDomain,D,R>
  {
    typedef GenericLocalFiniteElement< OrthonormalBasisFactory< dimDomain, SF, CF >,
        DGLocalCoefficientsFactory< OrthonormalBasisFactory< dimDomain, SF, CF > >,
        LocalL2InterpolationFactory< OrthonormalBasisFactory< dimDomain, SF, CF >,true >,
        dimDomain,D,R> Base;
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
