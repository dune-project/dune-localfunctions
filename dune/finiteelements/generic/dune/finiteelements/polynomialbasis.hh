// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_POLYNOMIALBASIS_HH
#define DUNE_POLYNOMIALBASIS_HH
#include <fstream>
#include <dune/common/fmatrix.hh>
#include <dune/finiteelements/coeffmatrix.hh>
#include <dune/finiteelements/monomialbasis.hh>
#include <dune/finiteelements/multiindex.hh>
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
  template<  class CM, class Eval >
  class PolynomialBasis
  {
    typedef PolynomialBasis< CM, Eval > This;

    typedef CM CoefficientMatrix;
    typedef Eval Evaluator;

    static const int dimension = Evaluator::dimension;
    typedef typename CoefficientMatrix::Field StorageField;

  public:
    typedef typename Evaluator::Basis Basis;
    typedef typename Evaluator::DomainVector DomainVector;

    PolynomialBasis (const Basis &basis,
                     const CoefficientMatrix &coeffMatrix,
                     int order, unsigned int size)
      : coeffMatrix_(&coeffMatrix),
        eval_(basis,order),
        size_(size)
    {
      assert(size <= coeffMatrix.size());
    }

    const unsigned int size () const
    {
      return size_;
    }

    template< class RangeVector >
    void evaluate ( const DomainVector &x,
                    std::vector< RangeVector > &values ) const
    {
      assert(values.size()>=size());
      /*
         eval_.evaluate( x );
         coeffMatrix_->mult( eval_, values );
       */
      eval_.evaluate( x );
      coeffMatrix_->mult( eval_.evaluate( x ), values );
    }

    template< class DomainVector, class RangeVector >
    void evaluate ( const DomainVector &x,
                    std::vector< RangeVector > &values ) const
    {
      typename This::DomainVector bx;
      for( int d = 0; d < dimension; ++d )
        field_cast( x[ d ], bx[ d ] );
      evaluate( bx, values );
    }

    template <class Topology,class FullMatrix>
    void printBasis(const std::string &name,
                    const FullMatrix &matrix)
    {
#if 0
      {
        std::ofstream out((name+".out").c_str());
        out.precision(15);
        out.setf(std::ios::scientific,std::ios::floatfield);
        matrix.print(out,size_);
      }
      {
        std::ofstream out((name+".gnu").c_str());
        out.precision(15);
        out.setf(std::ios::scientific,std::ios::floatfield);

        typedef Dune::MultiIndex<dimension> MI;
        typedef Dune::MonomialBasis< Topology, MI > MBasis;
        MBasis basis;
        const unsigned int size = basis.size( order_ );
        std::vector< MI > y( size );
        Dune::FieldVector< MI, dimension > x;
        for (int d=0; d<dimension; ++d)
          x[d].set(d);
        basis.evaluate( order_, x, &(y[0]) );
        coeffMatrix_->print(out,y,size_);
      }
#endif
    }
  protected:
    PolynomialBasis(const PolynomialBasis &);
    PolynomialBasis &operator=(const PolynomialBasis&);
    const CoefficientMatrix* coeffMatrix_;
    mutable Evaluator eval_;
    unsigned int size_;
  };

  /**
   * Specialized version of PolynomialBasis with FieldMatrix for matrix
   * coefficience and std::vector for container type with FieldVector as
   * value type. This class stores the coefficient matrix with can be
   * constructed via the fill method
   */
  template< class Eval >
  class PolynomialBasisWithMatrix
    : public PolynomialBasis<CoeffMatrix<typename Eval::Field> , Eval >
  {
    typedef typename Eval::Field StorageField;
    typedef Eval Evaluator;
    typedef CoeffMatrix< StorageField > CoefficientMatrix;

    typedef PolynomialBasisWithMatrix< Evaluator > This;
    typedef PolynomialBasis<CoefficientMatrix,Evaluator> Base;

    typedef typename Base::Basis Basis;
  public:
    PolynomialBasisWithMatrix (const Basis &basis,
                               int order)
      : Base(basis,coeffMatrix_,order,0)
    {}

    template <class FullMatrix>
    void fill(const FullMatrix& matrix)
    {
      coeffMatrix_.fill(matrix);
      this->size_ = matrix.rowSize();
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
