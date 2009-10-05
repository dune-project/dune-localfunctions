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

  template< int dimRange, class B, class SF >
  class PolynomialBasis
  {
    typedef PolynomialBasis< dimRange, B, SF > This;

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

    template< class RangeVector >
    void evaluate ( const DomainVector &x,
                    std::vector< RangeVector > &values ) const
    {
      basis_->evaluate( order_, x, basisEval_ );
      coeffMatrix_.mult( basisEval_, values );
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



  template< int dim, class SF,
      class Creator, class MBasis = VirtualMonomialBasis<dim,SF> >
  struct PolynomialBasisProvider
  {
    static const int dimension = dim;
    typedef SF Field;
    typedef PolynomialBasis<1,MBasis,SF> Basis;
    static const Basis &basis(unsigned int id,unsigned int order)
    {
      return instance().getBasis(id,order);
    }
  private:
    friend struct MonomialBasisProvider<dimension,SF>;
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
        MonomialBasisProvider<dimension,SF>::template callback<Creator>(id,order,basis_[order][id]);
      }
      else if (basis_[order][id] == 0)
        MonomialBasisProvider<dimension,SF>::template callback<Creator>(id,order,basis_[order][id]);
      return *(basis_[order][id]);
    }
    std::vector<FieldVector<Basis*,numTopologies> > basis_;
  };
}
#endif // DUNE_POLYNOMIALBASIS_HH
