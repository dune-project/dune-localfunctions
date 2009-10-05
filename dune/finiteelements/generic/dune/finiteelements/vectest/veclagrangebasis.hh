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
#include <dune/finiteelements/polynomialbasis.hh>
namespace Dune
{
  template <class Topology,class scalar_t,int dimR>
  struct VecLagrangeMatrix {
    static const unsigned int dimension = Topology::dimension;
    typedef Dune::AlgLib::Matrix< scalar_t > mat_t;
    typedef Dune::LagrangePointsCreator< scalar_t, dimension > LagrangePointsCreator;
    typedef LocalLagrangeInterpolationCreator< LagrangePointsCreator > LocalInterpolationCreator;
    typedef typename LocalInterpolationCreator::LocalInterpolation LocalInterpolation;

    explicit VecLagrangeMatrix( const unsigned int order )
    {
      Dune::MonomialBasis< Topology, scalar_t > basis( order );
      const LocalInterpolation &localInterpolation
        = LocalInterpolationCreator::template localInterpolation< Topology >( order );
      localInterpolation.interpolate( basis, matrix_ );
      LocalInterpolationCreator::release( localInterpolation );
      matrix_.invert();
    }
    unsigned int colSize(int row) const {
      return matrix_.cols()*dimR;
    }
    unsigned int rowSize() const {
      return matrix_.rows()*dimR;
    }
    const scalar_t operator() ( int r, int c ) const
    {
      assert(0<=r && r<rowSize());
      assert(0<=c && c<colSize(r));
      if (c%dimR == r%dimR)
        return matrix_(c/dimR,r/dimR);
      else
        return Zero<scalar_t>();
    }
    void print(std::ostream& out,int N = rowSize()) const {}
    mat_t matrix_;
  };
  template <class Topology,class scalar_t,int dimR>
  struct VecLagrangeMatrix1 {
    static const unsigned int dimension = Topology::dimension;
    typedef Dune::AlgLib::Matrix< scalar_t > mat_t;
    typedef Dune::LagrangePointsCreator< scalar_t, dimension > LagrangePointsCreator;
    typedef LocalLagrangeInterpolationCreator< LagrangePointsCreator > LocalInterpolationCreator;
    typedef typename LocalInterpolationCreator::LocalInterpolation LocalInterpolation;

    explicit VecLagrangeMatrix1( const unsigned int order )
      : matrix_(order)
    {
      baseSize = matrix_.rowSize()/dimR;
      mat_ = new double*[baseSize*dimR*dimR];
      int row = 0;
      for (int i=0; i<baseSize; ++i)
      {
        for (int r=0; r<dimR; ++r)
        {
          for (int rr=0; rr<dimR; ++rr)
          {
            mat_[row] = new double[baseSize];
            for (int j=0; j<baseSize; ++j)
            {
              mat_[row][j] = 0;
            }
            if (r==rr)
              mat_[row][i] = 1.;
            ++row;
          }
        }
      }
    }
    unsigned int colSize(int row) const {
      return baseSize;
    }
    unsigned int rowSize() const {
      return baseSize*dimR*dimR;
    }
    const scalar_t operator() ( int r, int c ) const
    {
      int rmod = r%dimR;
      int rdiv = r/dimR;
      scalar_t ret = 0;
      for (int k=0; k<baseSize*dimR; ++k) {
        assert(k*dimR+rmod<baseSize*dimR*dimR);
        assert(c<baseSize);
        ret += matrix_(rdiv,k)*mat_[k*dimR+rmod][c];
      }
      return ret;
    }
    VecLagrangeMatrix<Topology,scalar_t,dimR> matrix_;
    double **mat_;
    unsigned int baseSize;
  };

  // **************************************************

  template< unsigned int dim, int dimR, class SF, class CF = typename ComputeField< SF, 512 >::Type >
  struct VecLagrangeBasisCreator
  {
    static const unsigned int dimension = dim;

    typedef SF StorageField;

    typedef StorageField BasisField;

    typedef Dune::MonomialBasisProvider< dimension, BasisField > MonomialBasisProvider;
    typedef typename MonomialBasisProvider::Basis MonomialBasis;

    typedef AlgLib::MultiPrecision< Precision< CF >::value > ComputeField;
    typedef unsigned int Key;

    typedef VectorialEvaluator< MonomialBasis,dimR,value > Evaluator;
    typedef PolynomialBasisWithMatrix< Evaluator, SparseCoeffMatrix< StorageField, 1 > > Basis;
    // typedef VectorialEvaluator< MonomialBasis,1,value > Evaluator;
    // typedef PolynomialBasisWithMatrix< Evaluator, SparseCoeffMatrix< StorageField, dimR > > Basis;

    template< class Topology >
    static const Basis &basis ( const Key &order )
    {
      const MonomialBasis &monomialBasis = MonomialBasisProvider::template basis< Topology >( order );
      Basis *basis = new Basis( monomialBasis );

      VecLagrangeMatrix< Topology, ComputeField,dimR > matrix( order );
      // VecLagrangeMatrix1< Topology, ComputeField,dimR > matrix( order );

      basis->fill( matrix );
      return *basis;
    }

    static void release ( const Basis &basis )
    {
      delete &basis;
    }
  };
  #if 0
  template< int dim, int dimR, class SF, class CF >
  struct VecLagrangeBasisCreator
  {
    typedef VirtualMonomialBasis<dim,SF> MBasis;
    typedef SF StorageField;
    typedef AlgLib::MultiPrecision< Precision<CF>::value > ComputeField;
    static const int dimension = dim;
    typedef VectorialEvaluator<MBasis,dimR,derivative> Evaluator;
    typedef PolynomialBasisWithMatrix<Evaluator,SparseCoeffMatrix<StorageField> > Basis;
    typedef unsigned int Key;

    template <class Topology>
    struct Maker
    {
      static void apply(unsigned int order,Basis* &basis)
      {
        const MBasis &virtBasis = MonomialBasisProvider<dimension,StorageField>::template basis<Topology>(order);
        basis = new Basis(virtBasis);
        VecLagrangeMatrix<Topology,ComputeField,dimR> matrix(order);
        basis->fill(matrix);
        std::stringstream name;
        name << "lagrange_" << Topology::name() << "_p" << order;
        basis->template printBasis<Topology>(name.str(),matrix);
      }
    };
  };
  #endif
  template< int dim, int dimR, class SF, class CF = typename ComputeField< SF, 512 >::Type >
  struct VecLagrangeBasisProvider
    : public BasisProvider<VecLagrangeBasisCreator<dim,dimR,SF,CF> >
  {};
}
#endif // DUNE_ORTHONORMALBASIS_HH
