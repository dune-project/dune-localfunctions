// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_POLYNOMIALBASIS_HH
#define DUNE_POLYNOMIALBASIS_HH
#include <fstream>
#include <dune/finiteelements/coeffmatrix.hh>
#include <dune/finiteelements/monomialbasis.hh>
#include <dune/finiteelements/multiindex.hh>
namespace Dune
{

  // PolynomialBasis
  // ---------------

  template<  class B, class CM >
  class PolynomialBasis
  {
    typedef PolynomialBasis< B, CM > This;

    typedef B Basis;
    enum {dimension = Basis::dimension};

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
        typedef Dune::MonomialBasis< Topology, MI > Basis;
        Basis basis;
        const unsigned int size = basis.size( order_ );
        std::vector< MI > y( size );
        Dune::FieldVector< MI, dimension > x;
        for (int d=0; d<dimension; ++d)
          x[d].set(d);
        basis.evaluate( order_, x, &(y[0]) );
        coeffMatrix_->print(out,y,size_);
      }
    }
  protected:
    PolynomialBasis(const PolynomialBasis &);
    PolynomialBasis &operator=(const PolynomialBasis&);
    const Basis *basis_;
    const CoefficientMatrix* coeffMatrix_;
    mutable std::vector<StorageField> basisEval_;
    unsigned int order_,size_;
  };

  template< class B, class SF, int dimR >
  class PolynomialBasisWithMatrix
    : public PolynomialBasis<B,CoeffMatrix<FieldVector<SF,dimR> > >
  {
    typedef PolynomialBasisWithMatrix< B, SF, dimR > This;

    typedef B Basis;
    enum {dimension = Basis::dimension};

    static const int dimRange = dimR;
    typedef SF StorageField;
    typedef FieldVector<StorageField,dimRange> CoeffRangeVector;
    typedef CoeffMatrix< CoeffRangeVector > CoefficientMatrix;

    typedef PolynomialBasis<B,CoefficientMatrix> Base;

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


  template< class Creator >
  struct PolynomialBasisProvider
  {
    typedef typename Creator :: StorageField StorageField;
    static const int dimension = Creator :: dimension;
    typedef typename Creator :: Basis Basis;

    static const Basis &basis(unsigned int id,unsigned int order)
    {
      return instance().getBasis(id,order);
    }
  private:
    friend struct MonomialBasisProvider<dimension,StorageField>;
    enum { numTopologies = (1 << dimension) };
    PolynomialBasisProvider()
    {}
    static PolynomialBasisProvider &instance()
    {
      static PolynomialBasisProvider instance;
      return instance;
    }
    const Basis &getBasis(unsigned int id,unsigned int order)
    {
      if (order>=basis_.size())
      {
        basis_.resize(order+1,FieldVector<Basis*,numTopologies>(0));
        MonomialBasisProvider<dimension,StorageField>::template callback<Creator>(id,order,basis_[order][id]);
      }
      else if (basis_[order][id] == 0)
        MonomialBasisProvider<dimension,StorageField>::template callback<Creator>(id,order,basis_[order][id]);
      return *(basis_[order][id]);
    }
    std::vector<FieldVector<Basis*,numTopologies> > basis_;
  };
}
#endif // DUNE_POLYNOMIALBASIS_HH
