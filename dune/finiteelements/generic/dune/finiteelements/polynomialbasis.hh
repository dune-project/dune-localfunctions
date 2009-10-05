// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_MONOMIALBASIS_HH
#define DUNE_MONOMIALBASIS_HH
#include <vector>
namespace Dune
{
  template <class B,class M>
  class PolynomialBasis
  {
    typedef B Basis;
    typedef M Matrix;
    typedef PolynomialBasis<Basis,Matrix> This;

    Basis basis_;
    const Matrix& matrix_;
    std::vector<double> basisEval_;

  public:
    typedef typename Basis::Field Field;

    typedef typename Basis::DomainVector DomainVector;
    typedef typename Basis::RangeVector RangeVector;

    PolynomialBasis (const Matrix& matrix)
      : basis_(), matrix_(matrix), basisEval_(0)
    {}

    const unsigned int *sizes ( unsigned int order ) const
    {
      return basis_.sizes();
    }

    void evaluate ( const unsigned int order,
                    const DomainVector &x,
                    RangeVector *const values ) const
    {
      int size = basis_.sizes(order)[order];
      if (basisEval_.size() < size)
        basisEval_.resize(size);
      basis_.evaluate(order,x,&(basisEval_[0]));
      int posVec = 0;
      typename Matrix::Iterator rowEnd = matrix_.end();
      for ( typename Matrix::Iterator rowIt = matrix_.begin();
            rowIt != rowEnd; ++rowIt, ++posVec ) {
        values[ posVec ] = 0;
        int posMat = 0;
        typename Matrix::Row row = *rowIt;
        typename Matrix::Row::Iterator colEnd = row.end();
        for ( typename Matrix::Row::Iterator colIt = row.begin();
              colIt != colEnd; ++colIt, ++posMat ) {
          values[ posVec ] += (*colIt)*basisEval_[posMat];
        }
      }
    }

    void evaluate ( const unsigned int order,
                    const DomainVector &x,
                    std::vector< RangeVector > &values ) const
    {
      evaluate( order, x, &(values[ 0 ]) );
    }
  };
}
#endif // DUNE_MONOMIALBASIS_HH
