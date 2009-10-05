// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_POLYNOMIALBASIS_HH
#define DUNE_POLYNOMIALBASIS_HH
#include <fstream>
#include <dune/common/fmatrix.hh>
#include <dune/finiteelements/coeffmatrix.hh>
#include <dune/finiteelements/monomialbasis.hh>
#include <dune/finiteelements/multiindex.hh>
#include <dune/finiteelements/basisevaluator.hh>
namespace Dune
{

  // PolynomialBasis
  // ---------------

  /**
   * \tparam B Basis set with
   *           static const int dimension  -> dimension of reference element
   *           typedef DomainVector        -> coordinates in reference element
   *           int size(int order) const   -> number of basis functions
   *           void evaluate( order, x, val ) const
   *              int order
   *              DomainVector x
   *              Container val
   * \tparam CM stroage for coefficience with
   *           typedef Field -> field of coefficience
   *           static const int dimRange -> coeficience are of type
   *                                        FieldMatrix<Field,dimRange,dimRange>
   *           void mult( val, y )
   *              Container val
   *              std::vector<RangeVector> y
   * \tparam Container access to basis functions through forward iterator
   *           typedef value_type
   *           typedef const_iterator
   *           const_iterator begin()
   **/
  template< class Eval, class CM >
  class PolynomialBasis
  {
    typedef PolynomialBasis< Eval, CM > This;

    typedef CM CoefficientMatrix;
    typedef Eval Evaluator;
    typedef typename CoefficientMatrix::Field StorageField;
  public:
    static const int dimension = Evaluator::dimension;
    static const int dimRange = Evaluator::dimRange;
    typedef typename Evaluator::Basis Basis;
    typedef typename Evaluator::DomainVector DomainVector;
    typedef typename CoefficientMatrix::Field Field;

    PolynomialBasis (const Basis &basis,
                     const CoefficientMatrix &coeffMatrix,
                     unsigned int size)
      : coeffMatrix_(&coeffMatrix),
        eval_(basis),
        order_(basis.order()),
        size_(size)
    {
      assert(size_ <= coeffMatrix.size());
    }

    const unsigned int order () const
    {
      return order_;
    }
    const unsigned int size () const
    {
      return size_;
    }

    template< unsigned int deriv, class Vector >
    void evaluate ( const DomainVector &x,
                    Vector &values ) const
    {
      assert(values.size()>=size());
      coeffMatrix_->mult( eval_.template evaluate<deriv>( x ), values );
    }
    template< class Vector >
    void evaluate ( const DomainVector &x,
                    Vector &values ) const
    {
      evaluate<0>(x,values);
    }

    template< class DVector, class RVector >
    void evaluate ( const DVector &x,
                    RVector &values ) const
    {
      DomainVector bx;
      for( int d = 0; d < dimension; ++d )
        field_cast( x[ d ], bx[ d ] );
      evaluate<0>( bx, values );
    }

  protected:
    PolynomialBasis(const PolynomialBasis &);
    PolynomialBasis &operator=(const PolynomialBasis&);
    const CoefficientMatrix* coeffMatrix_;
    mutable Evaluator eval_;
    unsigned int order_,size_;
  };

  /**
   * Specialized version of PolynomialBasis with FieldMatrix for matrix
   * coefficience and std::vector for container type with FieldVector as
   * value type. This class stores the coefficient matrix with can be
   * constructed via the fill method
   */
  template< class Eval, class CM = SparseCoeffMatrix<typename Eval::Field,Eval::dimRange> >
  class PolynomialBasisWithMatrix
    : public PolynomialBasis< Eval, CM >
  {
    typedef typename Eval::Field StorageField;
    typedef Eval Evaluator;
    typedef CM CoefficientMatrix;

    typedef PolynomialBasisWithMatrix< Evaluator, CM > This;
    typedef PolynomialBasis<Evaluator,CM> Base;

    typedef typename Base::Basis Basis;
  public:
    PolynomialBasisWithMatrix (const Basis &basis)
      : Base(basis,coeffMatrix_,0)
    {}

    template <class FullMatrix>
    void fill(const FullMatrix& matrix)
    {
      coeffMatrix_.fill(matrix);
      this->size_ = coeffMatrix_.size();
    }
    template <class FullMatrix>
    void fill(const FullMatrix& matrix,int size)
    {
      coeffMatrix_.fill(matrix);
      assert(size<=coeffMatrix_.size());
      this->size_ = size;
    }

  private:
    PolynomialBasisWithMatrix(const PolynomialBasisWithMatrix &);
    PolynomialBasisWithMatrix &operator=(const PolynomialBasisWithMatrix &);
    CoefficientMatrix coeffMatrix_;
  };
}
#endif // DUNE_POLYNOMIALBASIS_HH
