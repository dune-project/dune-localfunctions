// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef BASISPRINT
#define BASISPRINT
#include <dune/localfunctions/utility/multiindex.hh>
#include <dune/localfunctions/utility/polynomialbasis.hh>
namespace Dune {
  /**********************************************
  * Methods for printing a PolynomialBasis.
  * Is achieved by using the MultiIndex class as
  * Field type and printing the results.
  * The basis and higher order derivatives can be
  * printed. This could be the basis for printing
  * routings providing C++ or matlab methods
  * for computing the basisfunctions for given
  * orders or reference elements.
  **********************************************/
  // default argument does not work for gcc 4.1.2
  // template <int deriv,class BasisFactory,class PrintField=typename BasisFactory::StorageField>
  template <int deriv,class BasisFactory,class PrintField,GeometryType::Id geometryId>
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

    MIBasis *miBasis = MIBasisFactory::template create<geometryId>( basis.basis().order());
    PrintBasis printBasis(*miBasis,basis.matrix(),basis.size());

    unsigned int size = printBasis.size();

    out << "% Number of base functions:  " << size << std::endl;
    out << "% Derivative order: " << deriv << std::endl;

    std::vector< FieldVector<
            FieldVector<Field,LFETensor<Field,dimension,deriv>::size>,
            PrintBasis::dimRange> > y( size );

    FieldVector< Field, dimension > x;
    for( int i = 0; i < dimension; ++i )
      x[ i ].set( i, 1 );
    printBasis.template evaluateSingle<deriv>( x, y );
    for (unsigned int i=0; i<size; ++i)
    {
      out << "$\\varphi_" << i << "(a,b,c)$&$=$&$" << std::endl;
      out << "( ";
      for (unsigned int r=0; r<PrintBasis::dimRange; ++r)
        out << y[i][r] << (r<PrintBasis::dimRange-1 ? " , $ \\\\ && $" : " )$ \\\\");
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
