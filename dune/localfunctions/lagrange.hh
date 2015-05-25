// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LAGRANGEFINITEELEMENT_HH
#define DUNE_LAGRANGEFINITEELEMENT_HH

#include <dune/localfunctions/utility/localfiniteelement.hh>
#include <dune/localfunctions/utility/dglocalcoefficients.hh>

#include <dune/localfunctions/lagrange/lagrangecoefficients.hh>
#include <dune/localfunctions/lagrange/interpolation.hh>
#include <dune/localfunctions/lagrange/lagrangebasis.hh>

namespace Dune
{
  /**
   * @brief Lagrange local finite elements for a given set of interpolation
   *        points.
   *
   * The class LP provides the points for the interpolation.
   * It has two template arguments, the first is the Field type to
   * use for evaluating the points the second the dimension
   * of the reference elements on which to construct the points.
   * It is instantiated with the desired order and has a template
   * method build taking a Topology to construct the points
   * (a std::vector of FieldVectors).
   * It also provides a static template method supports to indicate
   * if the point set can be build for a specified Topology.
   *
   * Examples include:
   * - EqualdistantPointSet: standard point set for lagrange points
   * - LobattoPointSet:      an approximate Freget type point set
   *                         (provided for simplex and generalized prism
   *                         topologies (i.e. not for a 3d pyramid)
   *
   * \ingroup Lagrange
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
  class LagrangeLocalFiniteElement
    : public GenericLocalFiniteElement< LagrangeBasisFactory< LP, dimDomain, SF, CF >,
          LagrangeCoefficientsFactory<LP, dimDomain, SF >,
          LagrangeInterpolationFactory< LP, dimDomain, SF > >
  {
    typedef GenericLocalFiniteElement< LagrangeBasisFactory< LP, dimDomain, SF, CF >,
        LagrangeCoefficientsFactory<LP, dimDomain, SF >,
        LagrangeInterpolationFactory< LP, dimDomain, SF > > Base;
  public:
    typedef typename Base::Traits Traits;

    /** \todo Please doc me !
     */
    LagrangeLocalFiniteElement ( const GeometryType &gt, unsigned int order )
      : Base( gt, order )
    {}
  };
}

#endif // #ifndef DUNE_LAGRANGEFINITEELEMENT_HH
