// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_POLYNOMIALORDER_HH
#define DUNE_POLYNOMIALORDER_HH

#include <dune/grid/genericgeometry/topologytypes.hh>

namespace Dune
{

  template< class Topology >
  struct PolynomialOrder;

  template<>
  struct PolynomialOrder< GenericGeometry::Point >
  {
    static unsigned int polynomialOrder ( unsigned int order )
    {
      return 0;
    }
  };

  template< class BaseTopology >
  struct PolynomialOrder< GenericGeometry::Prism< BaseTopology > >
  {
    static unsigned int polynomialOrder ( unsigned int order )
    {
      return PolynomialOrder< BaseTopology >::polynomialOrder( order ) + order;
    }
  };

  template< class BaseTopology >
  struct PolynomialOrder< GenericGeometry::Pyramid< BaseTopology > >
  {
    static unsigned int polynomialOrder ( unsigned int order )
    {
      return std::max( order, PolynomialOrder< BaseTopology >::polynomialOrder( order ) );
    }
  };



  template< class Topology >
  inline static unsigned int polynomialOrder ( unsigned int order )
  {
    return PolynomialOrder< Topology >::polynomialOrder( order );
  }

}

#endif // #ifndef DUNE_POLYNOMIALORDER_HH
