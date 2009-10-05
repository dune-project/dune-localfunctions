// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ORTHONORMALBASIS_HH
#define DUNE_ORTHONORMALBASIS_HH
#include <sstream>
#include <dune/finiteelements/polynomialbasis.hh>
#include <dune/finiteelements/basisprovider.hh>
#include <dune/finiteelements/basisprint.hh>
#include "orthonormalcompute.hh"
namespace Dune
{
  template< int dim, class SF, class CF >
  struct ONBasisCreator
  {
    typedef VirtualMonomialBasis<dim,SF> MBasis;
    typedef SF StorageField;
    typedef AlgLib::MultiPrecision< Precision<CF>::value > ComputeField;
    static const int dimension = dim;
    typedef SparseCoeffMatrix< StorageField > CoefficientMatrix;
    typedef StandardEvaluator<MBasis> Evaluator;
    typedef PolynomialBasis<Evaluator,CoefficientMatrix> Basis;
    typedef unsigned int Key;
    typedef typename GenericGeometry::SimplexTopology< dim >::type SimplexTopology;

    template <class Topology>
    struct Maker
    {
      static void apply(unsigned int order,Basis* &basis)
      {
        const MBasis &_basis = MonomialBasisProvider<dimension,StorageField>::template basis<SimplexTopology>(order);
        static CoefficientMatrix _coeffs;
        if ( _coeffs.size() <= _basis.size() )
        {
          ONB::ONBMatrix<Topology,ComputeField> matrix(order);
          _coeffs.fill(matrix);
          basis = new Basis(_basis,_coeffs,_basis.size());
        }
        else
          basis = new Basis(_basis,_coeffs,_basis.size());
        {
          typedef MultiIndex< dimension > MIField;
          typedef VirtualMonomialBasis<dim,MIField> MBasisMI;
          typedef PolynomialBasis<StandardEvaluator<MBasisMI>,SparseCoeffMatrix<StorageField> > BasisMI;
          const MBasisMI &_mBasisMI = MonomialBasisProvider<dimension,MIField>::template basis<SimplexTopology>(order);
          BasisMI basisMI(_mBasisMI,_coeffs,_basis.size());
          std::stringstream name;
          name << "onb_" << Topology::name() << "_p" << order;
          std::ofstream out(name.str().c_str());
          basisPrint<0>(out,basisMI);
        }
      }
    };
  };

  template< int dim, class SF, class CF = typename ComputeField< SF, 512 >::Type >
  struct OrthonormalBasisProvider
    : public BasisProvider<ONBasisCreator<dim,SF,CF> >
  {};
}
#endif // DUNE_ORTHONORMALBASIS_HH
