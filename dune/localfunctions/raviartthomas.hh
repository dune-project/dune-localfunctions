// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_RAVIARTTHOMASFINITEELEMENT_HH
#define DUNE_RAVIARTTHOMASFINITEELEMENT_HH

#include <dune/localfunctions/utility/localfiniteelement.hh>
#include <dune/localfunctions/raviartthomas/raviartthomasbasis.hh>

namespace Dune
{
  /**
   * @brief Raviart-Thomas basis functions
   *
   * These basis functions are at the moment only available
   * for simplex geometry types.
   *
   * \tparam dimDomain dimension of reference elements
   * \tparam D domain for basis functions
   * \tparam R range for basis functions
   * \tparam SF storage field for basis matrix
   * \tparam CF compute field for basis matrix
   **/
  template< unsigned int dimDomain, class D, class R,
      class SF=R, class CF=SF >
  class RaviartThomasLocalFiniteElement
    : public GenericLocalFiniteElement< RaviartThomasBasisFactory< dimDomain, SF, CF >,
          RaviartThomasCoefficientsFactory< dimDomain >,
          RaviartThomasL2InterpolationFactory< dimDomain, SF > >
  {
    typedef GenericLocalFiniteElement< RaviartThomasBasisFactory< dimDomain, SF, CF >,
        RaviartThomasCoefficientsFactory< dimDomain >,
        RaviartThomasL2InterpolationFactory< dimDomain, SF > > Base;
  public:
    using typename Base::Traits;

    /** \todo Please doc me */
    RaviartThomasLocalFiniteElement ( const GeometryType &gt, unsigned int order )
      : Base( gt, order )
    {}
  };

} // namespace Dune

#endif // #ifndef DUNE_RAVIARTTHOMASFINITEELEMENT_HH
