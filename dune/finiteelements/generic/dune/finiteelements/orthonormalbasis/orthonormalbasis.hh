// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ORTHONORMALBASIS_HH
#define DUNE_ORTHONORMALBASIS_HH
#include <sstream>
#include <dune/finiteelements/polynomialbasis.hh>
#include <dune/finiteelements/basisprovider.hh>
#include "orthonormalcompute.hh"
namespace Dune
{
  template< int dim, class SF, class CF >
  struct ONBasisCreator
  {
    typedef StandardMonomialBasis<dim,SF> MBasis;
    typedef SF StorageField;
    typedef AlgLib::MultiPrecision< Precision<CF>::value > ComputeField;
    static const int dimension = dim;
    typedef SparseCoeffMatrix< StorageField > CoefficientMatrix;
    typedef MonomialEvaluator<MBasis> Evaluator;
    typedef PolynomialBasis<Evaluator,CoefficientMatrix> Basis;
    typedef unsigned int Key;

    template <class Topology>
    struct Maker
    {
      static void apply(unsigned int order,Basis* &basis)
      {
        static MBasis _basis;
        static CoefficientMatrix _coeffs;
        if ( _coeffs.size() <= _basis.size(order) )
        {
          ONB::ONBMatrix<Topology,ComputeField> matrix(order);
          _coeffs.fill(matrix);
          basis = new Basis(_basis,_coeffs,order,_basis.size(order));
          std::stringstream name;
          name << "onb_" << Topology::name() << "_p" << order;
          basis->template printBasis<Topology>(name.str(),matrix);
        }
        else
          basis = new Basis(_basis,_coeffs,order,_basis.size(order));
      }
    };
  };

  template< int dim, class SF, class CF = typename ComputeField< SF, 512 >::Type >
  struct OrthonormalBasisProvider
    : public BasisProvider<ONBasisCreator<dim,SF,CF> >
  {};
}
#endif // DUNE_ORTHONORMALBASIS_HH
