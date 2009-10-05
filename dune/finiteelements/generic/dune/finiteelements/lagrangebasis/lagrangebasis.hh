// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LAGRANGEBASIS_HH
#define DUNE_LAGRANGEBASIS_HH

#include <fstream>
#include <dune/alglib/multiprecision.hh>
#include <dune/alglib/matrix.hh>

#include <dune/finiteelements/lagrangebasis/lagrangepoints.hh>
#include <dune/finiteelements/lagrangebasis/interpolation.hh>
#include <dune/finiteelements/basisprovider.hh>
#include <dune/finiteelements/basisprint.hh>
#include <dune/finiteelements/polynomialbasis.hh>

namespace Dune
{

  template <class Topology,class scalar_t>
  struct LagrangeMatrix
  {
    enum {dim = Topology::dimension};
    typedef Dune::AlgLib::Matrix< scalar_t > mat_t;
    LagrangeMatrix(int order)
    {
      Dune::MonomialBasis< Topology, scalar_t > basis(order);
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
    const Dune::FieldMatrix< scalar_t, 1,1 > operator() ( int r, int c ) const
    {
      return matrix_(c,r);
    }
    void print(std::ostream& out,int N = rowSize()) const {
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

  template< int dim, class SF, class CF >
  struct LagrangeBasisCreator
  {
    typedef SF StorageField;

    typedef StorageField BasisField;

    typedef VirtualMonomialBasis<dim,BasisField> MBasis;
    typedef StandardEvaluator<MBasis> Evaluator;

    typedef AlgLib::MultiPrecision< Precision<CF>::value > ComputeField;
    static const int dimension = dim;
    typedef PolynomialBasisWithMatrix<Evaluator,SparseCoeffMatrix<StorageField> > Basis;
    typedef unsigned int Key;

    template <class Topology>
    static void basis(unsigned int order,Basis* &basis)
    {
      const MBasis &virtBasis = MonomialBasisProvider<dimension,BasisField>::template basis<Topology>(order);

      basis = new Basis(virtBasis);
      LagrangeMatrix<Topology,ComputeField> matrix(order);
      basis->fill(matrix);
      {
        typedef MultiIndex< dimension > MIField;
        typedef VirtualMonomialBasis<dim,MIField> MBasisMI;
        typedef PolynomialBasisWithMatrix<StandardEvaluator<MBasisMI>,SparseCoeffMatrix<StorageField> > BasisMI;
        const MBasisMI &_mBasisMI = MonomialBasisProvider<dimension,MIField>::template basis<Topology>(order);
        BasisMI basisMI(_mBasisMI);
        basisMI.fill(matrix);
        std::stringstream name;
        name << "lagrange_" << Topology::name() << "_p" << order;
        std::ofstream out(name.str().c_str());
        basisPrint<1>(out,basisMI);
      }
    }
  };



  template< int dim, class SF, class CF = typename ComputeField< SF, 512 >::Type >
  struct LagrangeBasisProvider
    : public BasisProvider<LagrangeBasisCreator<dim,SF,CF> >
  {};

}
#endif // DUNE_ORTHONORMALBASIS_HH
