// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_POLYNOMIALBASIS_HH
#define DUNE_POLYNOMIALBASIS_HH

#include <fstream>
#include <numeric>

#include <dune/common/fmatrix.hh>

#include <dune/localfunctions/common/localbasis.hh>

#include <dune/localfunctions/utility/coeffmatrix.hh>
#include <dune/localfunctions/utility/monomialbasis.hh>
#include <dune/localfunctions/utility/multiindex.hh>
#include <dune/localfunctions/utility/basisevaluator.hh>

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
  class PolynomialBasis
  {
    typedef PolynomialBasis< Eval, CM > This;
    typedef Eval Evaluator;

  public:
    typedef CM CoefficientMatrix;

    typedef typename CoefficientMatrix::Field StorageField;

    static const unsigned int dimension = Evaluator::dimension;
    static const unsigned int dimRange = Evaluator::dimRange*CoefficientMatrix::blockSize;
    typedef LocalBasisTraits<D,dimension,FieldVector<D,dimension>,
        R,dimRange,FieldVector<R,dimRange>,
        FieldMatrix<R,dimRange,dimension> > Traits;
    typedef typename Evaluator::Basis Basis;
    typedef typename Evaluator::DomainVector DomainVector;
    template <class Fy>
    using HessianFyType = FieldVector<FieldMatrix<Fy,dimension,dimension>,dimRange>;
    using HessianType = HessianFyType<R>;

    PolynomialBasis (const Basis &basis,
                     const CoefficientMatrix &coeffMatrix,
                     unsigned int size)
      : basis_(basis),
        coeffMatrix_(&coeffMatrix),
        eval_(basis),
        order_(basis.order()),
        size_(size)
    {
      // assert(coeffMatrix_);
      // assert(size_ <= coeffMatrix.size()); // !!!
    }

    const Basis &basis () const
    {
      return basis_;
    }

    const CoefficientMatrix &matrix () const
    {
      return *coeffMatrix_;
    }

    unsigned int order () const
    {
      return order_;
    }

    unsigned int size () const
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

    //! \brief Evaluate Jacobian of all shape functions
    void evaluateHessian (const typename Traits::DomainType& x,         // position
                          std::vector<HessianType>& out) const          // return value
    {
      out.resize(size());
      hessian(x,out);
    }

    //! \brief Evaluate partial derivatives of all shape functions
    void partial (const std::array<unsigned int, dimension>& order,
                  const typename Traits::DomainType& in,                   // position
                  std::vector<typename Traits::RangeType>& out) const      // return value
    {
      out.resize(size());
      auto totalOrder = std::accumulate(order.begin(), order.end(), 0);
      if (totalOrder == 0) {
        evaluateFunction(in, out);
      }
      else if (totalOrder == 1) {
        std::vector<typename Traits::JacobianType> jacs(out.size());
        unsigned int k;
        for (unsigned int i=0;i<order.size();++i)
          if (order[i]==1) k=i;
        evaluateJacobian(in, jacs);
        for (unsigned int i=0;i<out.size();++i)
          for (unsigned int r=0;r<Traits::RangeType::dimension;++r)
            out[i][r] = jacs[i][r][k];
      }
      else if (totalOrder == 2) {
        std::vector<HessianType> hesss(out.size());
        int k=-1,l=-1;
        for (unsigned int i=0;i<order.size();++i) {
          if (order[i] >= 1 && k == -1)
            k = i;
          else if (order[i]==1) l=i;
        }
        if (l==-1) l=k;
        evaluateHessian(in, hesss);
        for (unsigned int i=0;i<out.size();++i)
          for (unsigned int r=0;r<Traits::RangeType::dimension;++r)
            out[i][r] = hesss[i][r][k][l];
      }
      else {
        DUNE_THROW(NotImplemented, "Desired derivative order is not implemented");
      }
    }

    template< unsigned int deriv, class F >
    void evaluate ( const DomainVector &x, F *values ) const
    {
      coeffMatrix_->mult( eval_.template evaluate<deriv>( x ), size(), values);
    }
    template< unsigned int deriv, class DVector, class F >
    void evaluate ( const DVector &x, F *values ) const
    {
      assert( DVector::dimension == dimension);
      DomainVector bx;
      for( int d = 0; d < dimension; ++d )
        field_cast( x[ d ], bx[ d ] );
      evaluate<deriv>( bx, values );
    }

    template <bool dummy,class DVector>
    struct Convert
    {
      static DomainVector apply( const DVector &x )
      {
        assert( DVector::dimension == dimension);
        DomainVector bx;
        for( unsigned int d = 0; d < dimension; ++d )
          field_cast( x[ d ], bx[ d ] );
        return bx;
      }
    };
    template <bool dummy>
    struct Convert<dummy,DomainVector>
    {
      static const DomainVector &apply( const DomainVector &x )
      {
        return x;
      }
    };
    template< unsigned int deriv, class DVector, class RVector >
    void evaluate ( const DVector &x, RVector &values ) const
    {
      assert(values.size()>=size());
      const DomainVector &bx = Convert<true,DVector>::apply(x);
      coeffMatrix_->mult( eval_.template evaluate<deriv>( bx ), values );
    }

    template <class Fy>
    void evaluate ( const DomainVector &x, std::vector<FieldVector<Fy,dimRange> > &values ) const
    {
      evaluate<0>(x,values);
    }
    template< class DVector, class RVector >
    void evaluate ( const DVector &x, RVector &values ) const
    {
      assert( DVector::dimension == dimension);
      DomainVector bx;
      for( unsigned int d = 0; d < dimension; ++d )
        field_cast( x[ d ], bx[ d ] );
      evaluate<0>( bx, values );
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
    void jacobian ( const DomainVector &x,
        std::vector<FieldMatrix<Fy,dimRange,dimension> > &values ) const
    {
      assert(values.size()>=size());
      evaluateSingle<1>(x,reinterpret_cast<std::vector<FieldVector<Fy,dimRange*dimension> >&>(values));
    }
    template< class DVector, class RVector >
    void jacobian ( const DVector &x, RVector &values ) const
    {
      assert( DVector::dimension == dimension);
      DomainVector bx;
      for( unsigned int d = 0; d < dimension; ++d )
        field_cast( x[ d ], bx[ d ] );
      jacobian( bx, values );
    }
    template <class Fy>
    void hessian ( const DomainVector &x,
        std::vector<HessianFyType<Fy>> &values ) const
    {
      assert(values.size()>=size());
      // only upper part of hessians matrix is computed - so we have
      // y[0] = FV< FV<Fy,d*(d+1)/2>, dimRange>
      const unsigned int hsize = LFETensor<Fy,dimension,2>::size;
      std::vector< FieldVector< FieldVector<Fy,hsize>, dimRange> > y( size() );
      evaluateSingle<2>(x, y);
      unsigned int q = 0;
      for (unsigned int i = 0; i < size(); ++i)
        for (unsigned int r = 0; r < dimRange; ++r)
        {
          q = 0;
          // tensor-based things follow unintuitive index sceme
          // e.g. for dim = 3, the k-l index of y is 00,01,11,02,12,22, i.e. partial derivatives
          // are ordered: xx,xy,yy,xz,yz,zz

          // Fill values 'directionwise'
          for (unsigned int k = 0; k < dimension; ++k)
            for (unsigned int l = 0; l <= k; ++l)
            {

              values[i][r][k][l] = y[i][r][q];
              values[i][r][l][k] = y[i][r][q];
              assert(q < hsize);
              ++q;
            }
        }
      // evaluateSingle<2>(x,reinterpret_cast<std::vector<FieldVector<Fy,dimRange*dimension*dimension> >&>(values));
    }
    template< class DVector, class HVector >
    void hessian ( const DVector &x, HVector &values ) const
    {
      assert( DVector::dimension == dimension);
      DomainVector bx;
      for( unsigned int d = 0; d < dimension; ++d )
        field_cast( x[ d ], bx[ d ] );
      hessian( bx, values );
    }

    template <class Fy>
    void integrate ( std::vector<Fy> &values ) const
    {
      assert(values.size()>=size());
      coeffMatrix_->mult( eval_.template integrate(), values );
    }

  protected:
    PolynomialBasis(const PolynomialBasis &other)
      : basis_(other.basis_),
        coeffMatrix_(other.coeffMatrix_),
        eval_(basis_),
        order_(basis_.order()),
        size_(other.size_)
    {}
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
    typedef PolynomialBasis<Evaluator, CM, D, R> Base;

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
