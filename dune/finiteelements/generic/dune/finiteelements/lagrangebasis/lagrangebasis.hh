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


  template< int dim, class SF, class CF = typename ComputeField< SF, 512 >::Type >
  struct LagrangeBasisProvider
  {
    static const int dimension = dim;
    typedef SF Field;
    typedef VirtualMonomialBasis<dimension,SF> VirtualBasis;
    typedef PolynomialBasis<1,VirtualBasis,SF> Basis;
    static const Basis &basis(unsigned int id,unsigned int order)
    {
      return instance().getBasis(id,order);
    }
  private:
    friend struct MonomialBasisProvider<dimension,SF>;
    enum { numTopologies = (1 << dimension) };
    LagrangeBasisProvider()
    {}
    static LagrangeBasisProvider &instance()
    {
      static LagrangeBasisProvider instance;
      return instance;
    }
    const Basis &getBasis(unsigned int id,unsigned int order)
    {
      if (order>=basis_.size())
      {
        basis_.resize(order+1,FieldVector<Basis*,numTopologies>(0));
        MonomialBasisProvider<dimension,SF>::callback(*this,id,order);
      }
      else if (basis_[order][id] == 0)
        MonomialBasisProvider<dimension,SF>::callback(*this,id,order);
      return *(basis_[order][id]);
    }
    template <class Topology>
    void apply(const VirtualBasis &basis,unsigned int order)
    {
      const unsigned int id = Topology::id;
      Basis *newBasis = new Basis(basis,order);
      basis_[order][id] = newBasis;
      LagrangeMatrix<Topology,CF> matrix(order);
      basis_[order][id]->fill(matrix);
    }
    std::vector<FieldVector<Basis*,numTopologies> > basis_;
  };
}
#endif // DUNE_ORTHONORMALBASIS_HH
