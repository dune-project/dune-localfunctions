// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ORTHONORMALBASIS_HH
#define DUNE_ORTHONORMALBASIS_HH
#include <sstream>
#include <dune/finiteelements/polynomialbasis.hh>
#include "orthonormalcompute.hh"
namespace Dune
{
  template< class CF >
  struct ONBasisCreator
  {
    template <class Topology,class VirtualBasis,class Basis>
    static void apply(const VirtualBasis &virtBasis,unsigned int order,Basis* &basis)
    {
      static StandardMonomialBasis<Topology::dimension,typename VirtualBasis::Field> _basis;
      basis = new Basis(_basis,order);
      ONB::ONBMatrix<Topology,CF> matrix(order);
      basis->fill(matrix);
    }
  };

  template< int dim, class SF, class CF = typename ComputeField< SF, 512 >::Type >
  struct OrthonormalBasisProvider
    : public PolynomialBasisProvider<dim,SF,ONBasisCreator<CF>,
          StandardMonomialBasis<dim,SF> >
  {};
}
#endif // DUNE_ORTHONORMALBASIS_HH
