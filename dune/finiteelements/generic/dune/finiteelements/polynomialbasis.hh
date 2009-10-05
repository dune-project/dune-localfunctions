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
  template< int dimRange, class B, class SF>
  class PolynomialBasis
  {
    typedef B Basis;
    enum {dimension = Basis::dimension};

    typedef SF StorageField;

  public:
    typedef typename Basis::DomainVector DomainVector;
    typedef FieldVector<StorageField,dimRange> CoeffRangeVector;

    PolynomialBasis (const Basis &basis, int order)
      : basis_(&basis),
        basisEval_(basis.size(order)), order_(order)
    {}

    const int size () const
    {
      return basis_->size(order_);
    }

    template <class RangeVector>
    void evaluate ( const DomainVector &x,
                    std::vector< RangeVector > &values ) const
    {
      basis_->evaluate(order_,x,basisEval_);
      coeffMatrix_.mult(basisEval_,values);
    }

    template <class DomainVector,class RangeVector>
    void evaluate ( const DomainVector &x,
                    std::vector< RangeVector > &values ) const
    {
      DomainVector bx;
      for (int d=0; d<dimension; ++d)
        field_cast(x[d], bx[ d ]);
      evaluate(bx,values);
    }

    void print(std::ofstream &out) const {
      /*
         typedef typename Basis::Topology Topology;
         typedef Dune::MultiIndex<dimension> MI;
         typedef Dune::MonomialBasis< Topology, MI > Basis;
         Basis basis;
         const unsigned int size = basis.size( order_ );
         std::vector< MI > y( size );
         Dune::FieldVector< MI, dimension > x;
         for (int d=0;d<dimension;++d)
          x[d].set(d);
         basis.evaluate( order_, x, &(y[0]) );
         coeffMatrix_.print(out,y);
       */
    }

    template <class FullMatrix>
    void fill(const FullMatrix& matrix)
    {
      coeffMatrix_.fill(matrix);
      {
        std::ofstream out("coeffs.out");
        out.precision(15);
        out.setf(std::ios::scientific,std::ios::floatfield);
        matrix.print(out);
      }
      {
        std::ofstream out("coeffs.gnu");
        out.precision(15);
        out.setf(std::ios::scientific,std::ios::floatfield);
        print(out);
      }
    }
  private:
    PolynomialBasis(const PolynomialBasis &);
    PolynomialBasis &operator=(const PolynomialBasis&);
    const Basis *basis_;
    mutable std::vector<StorageField> basisEval_;
    CoeffMatrix< CoeffRangeVector > coeffMatrix_;
    unsigned int order_;
  };
}
#endif // DUNE_POLYNOMIALBASIS_HH
