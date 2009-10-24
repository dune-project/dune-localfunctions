// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef BASISPRINT
#define BASISPRINT
#include <dune/finiteelements/generic/multiindex.hh>
#include <dune/finiteelements/generic/polynomialbasis.hh>
namespace Dune {
  template <int deriv,class BasisFactory>
  void basisPrint(std::ostream &out,
                  typename BasisFactory::Object &basis)
  {
    typedef typename BasisFactory::Object Basis;
    const int dimension = Basis::dimension;

    typedef MultiIndex< dimension > Field;
    typedef typename BasisFactory::template EvaluationBasisFactory<dimension,Field>::Type
    MIBasisFactory;
    typedef typename MIBasisFactory::Object MIBasis;
    typedef typename Basis::CoefficientMatrix CMatrix;
    typedef PolynomialBasis<StandardEvaluator<MIBasis>, CMatrix > PrintBasis;

    MIBasis *miBasis = MIBasisFactory::create(basis.basis().topologyId(),basis.basis().order());
    PrintBasis printBasis(*miBasis,basis.matrix(),basis.size());

    unsigned int size = printBasis.size();

    out << "% Number of base functions:  " << size << std::endl;
    out << "% Derivative order: " << deriv << std::endl;

    std::vector< Dune::FieldVector<Field,PrintBasis::dimRange> > y( size );
    FieldVector< Field, dimension > x;
    for( int i = 0; i < dimension; ++i )
      x[ i ].set( i, 1 );
    printBasis.template evaluate<deriv>( x, y );
    out << y << std::endl;
    MIBasisFactory::release(miBasis);
  }
};

#if GLFEM_BASIS_PRINT
{
  typedef MultiIndex< dimension > MIField;
  typedef VirtualMonomialBasis<dim,MIField> MBasisMI;
  typedef PolynomialBasisWithMatrix<StandardEvaluator<MBasisMI>,SparseCoeffMatrix<StorageField,dimension> > BasisMI;
  const MBasisMI &_mBasisMI = MonomialBasisProvider<dimension,MIField>::template basis<Topology>(order+1);
  BasisMI basisMI(_mBasisMI);
  basisMI.fill(matrix);
  std::stringstream name;
  name << "orthonormal_" << Topology::name() << "_p" << order << ".basis";
  std::ofstream out(name.str().c_str());
  basisPrint<0>(out,basisMI);
}
#endif

#endif // BASISPRINT
