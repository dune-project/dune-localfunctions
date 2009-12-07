// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_POLYNOMIALBASIS_HH
#define DUNE_POLYNOMIALBASIS_HH

#include <fstream>

#include <dune/common/fmatrix.hh>

#include <dune/finiteelements/common/localbasis.hh>

#include <dune/finiteelements/generic/common/coeffmatrix.hh>
#include <dune/finiteelements/generic/common/monomialbasis.hh>
#include <dune/finiteelements/generic/common/multiindex.hh>
#include <dune/finiteelements/generic/common/basisevaluator.hh>

namespace Dune
{

  // PolynomialBasis
  // ---------------

  /**
   * This is the basis class for a ''polynomial''
   * basis, i.e., a basis consisting of linear
   * combiniations of a underlying second basis set.
   * Examples are standard polynomials where the
   * underlying basis is given by the MonomialBasis
   * class. The basis evaluation is given by the matrix
   * vector multiplication between the coefficient
   * matrix and the vector filled by evaluating the
   * underlying basis set.
   * This class is constructed using a reference of
   * the underlying basis and the coefficient matrix.
   * A specialization holding an instance
   * of the coefficient matrix is provided by the class
   * template< class Eval, class CM = SparseCoeffMatrix<typename Eval::Field,Eval::dimRange> >
   * class PolynomialBasisWithMatrix;
   *
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
  template< class Eval, class CM, class D=double, class R=double >
  class PolynomialBasis :
    public C1LocalBasisInterface<
        C1LocalBasisTraits<D,Eval::dimension,FieldVector<D,Eval::dimension>,
            R,Eval::dimRange*CM::blockSize,FieldVector<R,Eval::dimRange*CM::blockSize>,
            FieldVector<FieldVector<R,Eval::dimension>,Eval::dimRange*CM::blockSize> >,
        PolynomialBasis<Eval,CM,D,R >
        >
  {
    typedef PolynomialBasis< Eval, CM > This;
    typedef Eval Evaluator;

  public:
    typedef CM CoefficientMatrix;

    typedef typename CoefficientMatrix::Field StorageField;

    static const unsigned int dimension = Evaluator::dimension;
    static const unsigned int dimRange = Evaluator::dimRange*CoefficientMatrix::blockSize;
    typedef C1LocalBasisTraits<D,dimension,FieldVector<D,dimension>,
        R,dimRange,FieldVector<R,dimRange>,
        FieldVector<FieldVector<R,dimension>,dimRange> > Traits;
    typedef typename Evaluator::Basis Basis;
    typedef typename Evaluator::DomainVector DomainVector;

    PolynomialBasis (const Basis &basis,
                     const CoefficientMatrix &coeffMatrix,
                     unsigned int size)
      : basis_(basis),
        coeffMatrix_(&coeffMatrix),
        eval_(basis),
        order_(basis.order()),
        size_(size)
    {
      assert(size_ <= coeffMatrix.size());
    }

    const Basis &basis () const
    {
      return basis_;
    }

    const CoefficientMatrix &matrix () const
    {
      return *coeffMatrix_;
    }

    const unsigned int order () const
    {
      return order_;
    }

    const unsigned int size () const
    {
      return size_;
    }

    //! \brief Evaluate all shape functions
    void evaluateFunction (const typename Traits::DomainType& x,
                           std::vector<typename Traits::RangeType>& out) const
    {
      out.resize(size());
      evaluate(x,out);
    }

    //! \brief Evaluate Jacobian of all shape functions
    void evaluateJacobian (const typename Traits::DomainType& x,         // position
                           std::vector<typename Traits::JacobianType>& out) const      // return value
    {
      out.resize(size());
      jacobian(x,out);
    }

    template< unsigned int deriv, class F >
    void evaluate ( const DomainVector &x, F *values ) const
    {
      coeffMatrix_->mult( eval_.template evaluate<deriv>( x ),
                          size(), values);
    }
    template< unsigned int deriv, class Vector >
    void evaluate ( const DomainVector &x, Vector &values ) const
    {
      assert(values.size()>=size());
      coeffMatrix_->mult( eval_.template evaluate<deriv>( x ), values );
    }
    template< unsigned int deriv, class Vector >
    void evaluateSingle ( const DomainVector &x, Vector &values ) const
    {
      assert(values.size()>=size());
      coeffMatrix_->template mult<deriv>( eval_.template evaluate<deriv>( x ), values );
    }
    template< unsigned int deriv, class Fy >
    void evaluateSingle ( const DomainVector &x,
                          std::vector< FieldVector<FieldVector<Fy,LFETensor<Fy,dimension,deriv>::size>,dimRange> > &values) const
    {
      evaluateSingle<deriv>(x,reinterpret_cast<std::vector< FieldVector<Fy,LFETensor<Fy,dimension,deriv>::size*dimRange> >&>(values));
    }
    template< unsigned int deriv, class Fy >
    void evaluateSingle ( const DomainVector &x,
                          std::vector< FieldVector<LFETensor<Fy,dimension,deriv>,dimRange> > &values) const
    {
      evaluateSingle<deriv>(x,reinterpret_cast<std::vector< FieldVector<Fy,LFETensor<Fy,dimension,deriv>::size*dimRange> >&>(values));
    }
    template <class Fy>
    void jacobian ( const DomainVector &x, std::vector<FieldMatrix<Fy,dimRange,dimension> > &values ) const
    {
      assert(values.size()>=size());
      evaluateSingle<1>(x,reinterpret_cast<std::vector<FieldVector<Fy,dimRange*dimension> >&>(values));
    }
    template <class Fy>
    void evaluate ( const DomainVector &x, std::vector<FieldVector<Fy,dimRange> > &values ) const
    {
      evaluateSingle<0>(x,values);
    }

    template< class DVector, class RVector >
    void evaluate ( const DVector &x, RVector &values ) const
    {
      assert( DVector::size == dimension);
      DomainVector bx;
      for( unsigned int d = 0; d < dimension; ++d )
        field_cast( x[ d ], bx[ d ] );
      evaluate<0>( bx, values );
    }
    template< unsigned int deriv, class DVector, class F >
    void evaluate ( const DVector &x, F *values ) const
    {
      assert( DVector::size == dimension);
      DomainVector bx;
      for( int d = 0; d < dimension; ++d )
        field_cast( x[ d ], bx[ d ] );
      evaluate<deriv>( bx, values );
    }
    template< unsigned int deriv, class DVector, class RVector >
    void evaluate ( const DVector &x, RVector &values ) const
    {
      assert( DVector::size == dimension);
      DomainVector bx;
      for( int d = 0; d < dimension; ++d )
        field_cast( x[ d ], bx[ d ] );
      evaluate<deriv>( bx, values );
    }

    template <class Fy>
    void integrate ( std::vector<Fy> &values ) const
    {
      assert(values.size()>=size());
      coeffMatrix_->mult( eval_.template integrate(), values );
    }

  protected:
    PolynomialBasis(const PolynomialBasis &);
    PolynomialBasis &operator=(const PolynomialBasis&);
    const Basis &basis_;
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
  template< class Eval, class CM = SparseCoeffMatrix<typename Eval::Field,Eval::dimRange>,
      class D=double, class R=double>
  class PolynomialBasisWithMatrix
    : public PolynomialBasis< Eval, CM, D, R >
  {
  public:
    typedef CM CoefficientMatrix;

  private:
    typedef Eval Evaluator;

    typedef PolynomialBasisWithMatrix< Evaluator, CM > This;
    typedef PolynomialBasis<Evaluator,CM> Base;

  public:
    typedef typename Base::Basis Basis;

    PolynomialBasisWithMatrix (const Basis &basis)
      : Base(basis,coeffMatrix_,0)
    {}

    template <class Matrix>
    void fill(const Matrix& matrix)
    {
      coeffMatrix_.fill(matrix);
      this->size_ = coeffMatrix_.size();
    }
    template <class Matrix>
    void fill(const Matrix& matrix,int size)
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