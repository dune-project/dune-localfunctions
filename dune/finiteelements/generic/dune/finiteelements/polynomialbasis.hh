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
  template<  class B, class CM, class Container >
  class PolynomialBasis
  {
    typedef PolynomialBasis< B, CM, Container > This;

    typedef B Basis;
    static const int dimension = Basis::dimension;

    typedef CM CoefficientMatrix;
    typedef typename CoefficientMatrix::Field StorageField;
    static const int dimRange = CoefficientMatrix::dimension;

  public:
    typedef typename Basis::DomainVector DomainVector;

    PolynomialBasis (const Basis &basis,
                     const CoefficientMatrix &coeffMatrix,
                     int order,
                     int size = 0)
      : basis_(&basis),
        coeffMatrix_(&coeffMatrix),
        basisEval_(basis.size(order)),
        order_(order),
        size_(size)
    { }

    const int size () const
    {
      return size_;
    }

    template< class RangeVector >
    void evaluate ( const DomainVector &x,
                    std::vector< RangeVector > &values ) const
    {
      assert(values.size()>=size());
      basis_->evaluate( order_, x, basisEval_ );
      coeffMatrix_->mult( basisEval_, values );
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
    const Basis *basis_;
    const CoefficientMatrix* coeffMatrix_;
    mutable Container basisEval_;
    unsigned int order_,size_;
  };

  /**
   * Specialized version of PolynomialBasis with FieldMatrix for matrix
   * coefficience and std::vector for container type with FieldVector as
   * value type. This class stores the coefficient matrix with can be
   * constructed via the fill method
   */
  template< class B, class SF, int dimR, class Container =  std::vector<Dune::FieldVector<SF,dimR> > >
  class PolynomialBasisWithMatrix
    : public PolynomialBasis<B,CoeffMatrix<FieldMatrix<SF,dimR,dimR> > , Container>
  {
    typedef PolynomialBasisWithMatrix< B, SF, dimR, Container > This;

    typedef B Basis;
    enum {dimension = Basis::dimension};

    static const int dimRange = dimR;
    typedef SF StorageField;
    typedef FieldMatrix<StorageField,dimRange,dimRange> CoeffRangeVector;
    typedef CoeffMatrix< CoeffRangeVector > CoefficientMatrix;

    typedef PolynomialBasis<Basis,CoefficientMatrix,Container> Base;

  public:
    typedef typename Basis::DomainVector DomainVector;

    PolynomialBasisWithMatrix (const Basis &basis,
                               int order)
      : Base(basis,coeffMatrix_,order)
    {}

    template <class FullMatrix>
    void fill(const FullMatrix& matrix)
    {
      coeffMatrix_.fill(matrix);
      this->size_ = matrix.rowSize();
    }

  private:
    PolynomialBasisWithMatrix(const PolynomialBasisWithMatrix &);
    PolynomialBasisWithMatrix &operator=(const PolynomialBasisWithMatrix &);
    CoeffMatrix< CoeffRangeVector > coeffMatrix_;
  };
}
#endif // DUNE_POLYNOMIALBASIS_HH
