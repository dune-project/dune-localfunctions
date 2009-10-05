// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LAGRANGEBASIS_HH
#define DUNE_LAGRANGEBASIS_HH
#include <fstream>
#include <dune/alglib/multiprecision.hh>
#include <dune/alglib/matrix.hh>

#include <dune/finiteelements/lagrangepoints.hh>
#include <dune/finiteelements/lagrangeinterpolation.hh>
#include <dune/finiteelements/polynomialbasis.hh>
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
    const Dune::FieldVector< scalar_t, 1 > operator() ( int r, int c ) const
    {
      return matrix_(c,r);
    }
    void print(std::ostream& out) const {
      int N = rowSize();
      for (int i=0; i<N; ++i) {
        out << "Polynomial : " << i << std::endl;
        for (int j=0; j<colSize(i); j++) {
          double v = matrix_(j,i).toDouble();
          if (fabs(v)<1e-20)
            out << 0 << "\t\t" << std::flush;
          else {
            Dune::AlgLib::MultiPrecision<128> v = matrix_(j,i);
            out << v << "\t\t" << std::flush;
          }
        }
        out << std::endl;
      }
    }
    mat_t matrix_;
  };


  template< class Topology, class SF, class CF = typename ComputeField< SF, 512 >::Type >
  struct LagrangeBasis
    : public PolynomialBasis<1,MonomialBasis<Topology,SF>,SF,CF>
  {
    typedef PolynomialBasis<1,MonomialBasis<Topology,SF>,SF,CF> Base;
    LagrangeBasis (int order)
      : Base(order)
    {
      LagrangeMatrix<Topology,CF> matrix(order);
      this->fill(matrix);
    }
  };
}
#endif // DUNE_ORTHONORMALBASIS_HH
