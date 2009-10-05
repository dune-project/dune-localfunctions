// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LAGRANGEBASIS_HH
#define DUNE_LAGRANGEBASIS_HH
#include <fstream>
#include <dune/alglib/multiprecision.hh>
#include <dune/alglib/matrix.hh>

#include <dune/finiteelements/lagrangepoints.hh>
#include <dune/finiteelements/monomialbasis.hh>
#include <dune/finiteelements/lagrangeinterpolation.hh>
#include <dune/finiteelements/coeffmatrix.hh>
#include <dune/finiteelements/multiindex.hh>
namespace Dune
{
  template <class Topology,class scalar_t>
  struct LagrangeMatrix {
    enum {dim = Topology::dimension};
    typedef Dune::AlgLib::Matrix< scalar_t > mat_t;
    LagrangeMatrix(int order)
    {
      Dune::MonomialBasis< Topology, scalar_t > basis;
      Dune::LocalLagrangeInterpolation< Topology, scalar_t  > interpolation( order );
      interpolation.interpolate( basis, matrix_ );
      matrix_.invert();
    }
    int colSize(int row) const {
      return matrix_.cols();
    }
    int rowSize() const {
      return matrix_.rows();
    }
    void set(int r,int c,scalar_t &v) const {
      v = matrix_(c,r);
    }
    void set(int r,int c,double &v) const {
      v = matrix_(c,r).toDouble();
    }
    void set(int r,int c,std::string &v) const {
      v = amp::ampf<128>(matrix_(c,r)).toDec();
    }
    void print(std::ostream& out) {
      int N = rowSize();
      for (int i=0; i<N; ++i) {
        out << "Polynomial : " << i << std::endl;
        for (int j=0; j<colSize(i); j++) {
          double v = 0;
          set(i,j,v);
          if (fabs(v)<1e-20)
            out << 0 << "\t\t" << std::flush;
          else {
            std::string v;
            set(i,j,v);
            out << v << "\t\t" << std::flush;
          }
        }
        out << std::endl;
      }
    }
    mat_t matrix_;
  };


  template <class Topology,class F>
  class LagrangeBasis
  {
    enum {dim = Topology::dimension};
    typedef LagrangeBasis<Topology,F> This;
    typedef StandardMonomialBasis<dim,F> Basis;
    static const unsigned int Precision = 1024;
    typedef Dune::AlgLib::MultiPrecision< Precision > scalar_t;

  public:
    typedef F Field;

    typedef typename Basis::DomainVector DomainVector;
    typedef typename Basis::RangeVector RangeVector;

    LagrangeBasis (int order)
      : basis_(), order_(order),
        basisEval_(basis_.size(order))
    {
      LagrangeMatrix<Topology,scalar_t> matrix(order);
      coeffMatrix_.fill(matrix);
      std::ofstream out("coeffs.out");
      out.precision(15);
      out.setf(std::ios::scientific,std::ios::floatfield);
      matrix.print(out);
      out << " ************ " << std::endl;
      print(out);
    }

    const int size () const
    {
      return basis_.size(order_);
    }

    void evaluate ( const DomainVector &x,
                    std::vector< RangeVector > &values ) const
    {
      basis_.evaluate(order_,x,basisEval_);
      coeffMatrix_.mult(basisEval_,values);
    }

    void print(std::ofstream &out) {
      typedef Dune::MultiIndex<dim> MI;
      typedef Dune::MonomialBasis< Topology, MI  > Basis;
      Basis basis;
      const unsigned int size = basis.size( order_ );
      std::vector< Dune::FieldVector< MI,1> > y( size );
      Dune::FieldVector< MI, dim > x;
      for (int d=0; d<dim; ++d)
        x[d].set(d);
      basis.evaluate( order_, x, y );
      coeffMatrix_.print(out,y);
    }
  private:
    const Basis basis_;
    int order_;
    mutable std::vector<RangeVector> basisEval_;
    CoeffMatrix<scalar_t> coeffMatrix_;
  };
}
#endif // DUNE_ORTHONORMALBASIS_HH
