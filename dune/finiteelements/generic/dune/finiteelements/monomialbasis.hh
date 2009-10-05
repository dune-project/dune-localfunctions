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

}

#endif
