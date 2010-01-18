// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DGLAGRANGEFINITEELEMENT_HH
#define DUNE_DGLAGRANGEFINITEELEMENT_HH

#include <dune/localfunctions/utility/localfiniteelement.hh>
#include <dune/localfunctions/utility/dglocalcoefficients.hh>
#include <dune/localfunctions/generic/lagrangebasis/lagrangecoefficients.hh>
#include <dune/localfunctions/generic/lagrangebasis/interpolation.hh>
#include <dune/localfunctions/generic/lagrangebasis/lagrangebasis.hh>

namespace Dune
{
  /**
   * @brief a dg space using lagrange basis functions
   *
   * Similar to the Dune::LagrangeLocalFiniteElement class
   * but this class uses the Dune::DGLocalCoefficients factory
   * to build a discontinuous space from the lagrange basis functions.
   * The point set used to define the lagrange points is given
   * as first template argument; examples are
   * - Dune::EquidistantPointSet
   * - Dune::LobattoPointSet
   *
   * \tparam LP a template class defining the points for the lagrange interpolation
   * \tparam dimDomain dimension of reference elements
   * \tparam D domain for basis functions
   * \tparam R range for basis functions
   * \tparam SF storage field for basis matrix
   * \tparam CF compute field for basis matrix
   **/
  template< template <class,unsigned int> class LP,
      unsigned int dimDomain, class D, class R,
      class SF=R, class CF=SF >
  class DGLagrangeLocalFiniteElement
    : public GenericLocalFiniteElement< LagrangeBasisFactory< LP, dimDomain, SF, CF >,
          DGLocalCoefficientsFactory< LagrangeBasisFactory< LP, dimDomain, SF, CF > >,
          LagrangeInterpolationFactory< LP, dimDomain, SF >,
          dimDomain,D,R>
  {
    typedef GenericLocalFiniteElement< LagrangeBasisFactory< LP, dimDomain, SF, CF >,
        DGLocalCoefficientsFactory< LagrangeBasisFactory< LP, dimDomain, SF, CF > >,
        LagrangeInterpolationFactory< LP, dimDomain, SF >,
        dimDomain,D,R> Base;
  public:
    typedef typename Base::Traits Traits;

    /** \todo Please doc me !
     */
    DGLagrangeLocalFiniteElement ( unsigned int topologyId, unsigned int order )
      : Base( topologyId, order )
    {}
  };

}

#endif // #ifndef DUNE_DGLAGRANGEFINITEELEMENT_HH
