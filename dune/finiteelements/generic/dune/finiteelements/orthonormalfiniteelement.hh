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

  template< unsigned int dimDomain, class D, class R,
      class SF=D, class CF=D >

  class OrthonormalLocalFiniteElement
    : public GenericLocalFiniteElement< OrthonormalBasisProvider< dimDomain, SF, CF >,
          DGLocalCoefficientsCreator< OrthonormalBasisProvider< dimDomain, SF, CF > >,
          LocalL2InterpolationCreator< OrthonormalBasisProvider< dimDomain, SF, CF > >,
          dimDomain,D,R>
  {
    typedef GenericLocalFiniteElement< OrthonormalBasisProvider< dimDomain, SF, CF >,
        DGLocalCoefficientsCreator< OrthonormalBasisProvider< dimDomain, SF, CF > >,
        LocalL2InterpolationCreator< OrthonormalBasisProvider< dimDomain, SF, CF > >,
        dimDomain,D,R> Base;
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
