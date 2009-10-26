// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef BASISPRINT
#define BASISPRINT
#include <dune/finiteelements/generic/multiindex.hh>
#include <dune/finiteelements/generic/polynomialbasis.hh>
namespace Dune {
  template <int deriv,class BasisFactory,class PrintField=typename BasisFactory::StorageField>
  void basisPrint(std::ostream &out,
                  typename BasisFactory::Object &basis)
  {
    typedef typename BasisFactory::Object Basis;
    const int dimension = Basis::dimension;

    typedef MultiIndex< dimension, PrintField > Field;
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

    /*
        std::vector< FieldVector<
           LFETensor<Field,dimension,deriv>,PrintBasis::dimRange> >
          y( size );
     */
    std::vector< FieldVector<
            FieldVector<Field,LFETensor<Field,dimension,deriv>::size>,
            PrintBasis::dimRange> > y( size );

    FieldVector< Field, dimension > x;
    for( int i = 0; i < dimension; ++i )
      x[ i ].set( i, 1 );
    printBasis.template evaluateSingle<deriv>( x, y );
    for (unsigned int i=0; i<size; ++i)
    {
      out << "func_" << i << ":" << std::endl;
      out << "( ";
      for (unsigned int r=0; r<PrintBasis::dimRange; ++r)
        out << y[i][r] << (r<PrintBasis::dimRange-1 ? " , " : " )");
      out << std::endl;
    }
    MIBasisFactory::release(miBasis);
  }
  template <int deriv,class BasisFactory,class PrintField=typename BasisFactory::StorageField>
  void basisPrint(std::ostream &out,
                  typename BasisFactory::Key &key)
  {
    typename BasisFactory::Object *basis = BasisFactory::create(key);
    basisPrint<deriv,BasisFactory,PrintField>(out,*basis);
    BasisFactory::release(basis);
  }
}


#endif // BASISPRINT
