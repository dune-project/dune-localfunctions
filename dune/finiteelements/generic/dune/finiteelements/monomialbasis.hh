// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_MONOMIALBASIS_HH
#define DUNE_MONOMIALBASIS_HH

#include <dune/common/fvector.hh>

#include <dune/grid/genericgeometry/topologytypes.hh>

namespace Dune
{

  // MonomialBasis
  // -------------

  template< class Topology, class F, unsigned int k >
  class MonomialBasis;

  template< class F, unsigned int k >
  class MonomialBasis< GenericGeometry::Point, F, k >
  {
    typedef MonomialBasis< GenericGeometry::Point, F, k > This;

  public:
    typedef GenericGeometry::Point Topology;
    typedef F Field;

    static const unsigned int dimDomain = Topology::dimension;
    static const unsigned int dimRange = 1;

    typedef FieldVector< Field, dimDomain > DomainVector;
    typedef FieldVector< Field, dimRange > RangeVector;

    static const unsigned int order = k;
    static const unsigned int numBasisFunctions = 1;
  };

  template< class BaseTopology, class F, unsigned int k >
  class MonomialBasis< GenericGeometry::Prism< BaseTopology >, F, k >
  {
    typedef MonomialBasis< GenericGeometry::Prism< BaseTopology >, F, k > This;

  public:
    typedef GenericGeometry::Prism< BaseTopology > Topology;
    typedef F Field;

    static const unsigned int dimDomain = Topology::dimension;
    static const unsigned int dimRange = 1;

    typedef FieldVector< Field, dimDomain > DomainVector;
    typedef FieldVector< Field, dimRange > RangeVector;

  public:
    static const unsigned int order = k;

    // number of basis functions with exacly order k
    static const unsigned int numBasisFunctions
      = (order+1) * MonomialBasis< BaseTopology, Field, order >::numBasisFunctions
        + MonomialBasis< Topology, Field, order-1 >::numBasisFunctions;
  };

  template< class BaseTopology, class F, unsigned int k >
  class MonomialBasis< GenericGeometry::Pyramid< BaseTopology >, F, k >
  {
    typedef MonomialBasis< GenericGeometry::Pyramid< BaseTopology >, F, k > This;

  public:
    typedef GenericGeometry::Pyramid< BaseTopology > Topology;
    typedef F Field;

    static const unsigned int dimDomain = Topology::dimension;
    static const unsigned int dimRange = 1;

    typedef FieldVector< Field, dimDomain > DomainVector;
    typedef FieldVector< Field, dimRange > RangeVector;

  public:
    static const unsigned int order = k;

    // number of basis functions with exacly order k
    static const unsigned int numBasisFunctions
      = MonomialBasis< BaseTopology, Field, order >::numBasisFunctions
        + MonomialBasis< Topology, Field, order-1 >::numBasisFunctions;
  };
}

#endif
