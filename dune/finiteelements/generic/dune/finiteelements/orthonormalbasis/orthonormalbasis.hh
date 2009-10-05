// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ORTHONORMALBASIS_HH
#define DUNE_ORTHONORMALBASIS_HH
#include <sstream>
#include <dune/finiteelements/polynomialbasis.hh>
#include "orthonormalcompute.hh"
namespace Dune
{

  template< class Topology, class SF, class CF = typename ComputeField< SF, 512 >::Type >
  class OrthonormalBasis
    : public PolynomialBasis<1,StandardMonomialBasis<Topology::dimension,SF>,SF>
  {
    enum {dimension = Topology::dimension};

    typedef SF StorageField;
    typedef CF ComputationField;

    typedef StandardMonomialBasis<dimension,StorageField> Basis;
    typedef PolynomialBasis<1,Basis,SF> Base;

  public:
    typedef typename Basis::DomainVector DomainVector;

    OrthonormalBasis (int order)
      : Base(order)
    {
      ONB::ONBMatrix<Topology,ComputationField> onbMatrix(order);
      fill(onbMatrix);
    }
  };
}
#endif // DUNE_ORTHONORMALBASIS_HH
