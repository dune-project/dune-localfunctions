// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LAGRANGEBASIS_HH
#define DUNE_LAGRANGEBASIS_HH
#include <fstream>
#include <dune/alglib/multiprecision.hh>
#include <dune/alglib/matrix.hh>

#include <dune/finiteelements/lagrangepoints.hh>
#include <dune/finiteelements/lagrangeinterpolation.hh>
#include <dune/finiteelements/basisprovider.hh>
#include <dune/finiteelements/polynomialbasis.hh>
namespace Dune
{
  template <class Topology,class scalar_t,int dimR>
  struct VecLagrangeMatrix {
    enum {dim = Topology::dimension};
    typedef Dune::AlgLib::Matrix< scalar_t > mat_t;
    VecLagrangeMatrix(int order)
    {
      Dune::MonomialBasis< Topology, scalar_t > basis;
      Dune::LocalLagrangeInterpolation< Topology, scalar_t  > interpolation( order );
      interpolation.interpolate( basis, matrix_ );
      matrix_.invert();
    }
    int colSize(int row) const {
      return matrix_.cols()*dimR;
    }
    int rowSize() const {
      return matrix_.rows()*dimR;
    }
    const Dune::FieldMatrix< scalar_t, dimR,dimR > operator() ( int r, int c ) const
    {
      Dune::FieldMatrix< scalar_t, dimR,dimR > ret(0);
      ret[r%dimR][r%dimR] = matrix_(c/dimR,r/dimR);
      return ret;
    }
    void print(std::ostream& out,int N = rowSize()) const {}
    mat_t matrix_;
  };

  template <int d,class F>
  struct VecContainer
  {
    struct Iterator
    {
      Iterator(const std::vector<F>& mBasis, bool end)
        : mBasis_(mBasis),
          pos_(0), dim_(0),
          done_(end) {
        val_[dim_] = mBasis_[pos_];
      }
      Iterator(const Iterator& other)
        : mBasis_(other.mBasis_),
          pos_(other.pos_), dim_(other.dim_),
          done_(other.done_)
      {}
      Iterator &operator=(const Iterator& other)
      {
        mBasis_ = other.mBasis;
        pos_ = other.pos_;
        dim_ = other.dim_;
        done_ = other.done_;
      }
      const FieldVector<F,d>& operator*() const
      {
        assert(!done_);
        return val_;
      }
      const Iterator &operator++()
      {
        val_[dim_] = 0;
        if (dim_ == d-1)
        {
          dim_ = 0;
          ++pos_;
        }
        else
          ++dim_;
        if (pos_==mBasis_.size())
          done_ = true;
        else
          val_[dim_] = mBasis_[pos_];
        return *this;
      }
      FieldVector<F,d> val_;
      const std::vector<F>& mBasis_;
      unsigned int pos_;
      int dim_;
      bool done_;
    };
    typedef Dune::FieldVector<F,d> value_type;
    typedef Iterator const_iterator;
    VecContainer(int size)
      : mBasis_(size) {}
    Iterator begin() const
    {
      return Iterator(mBasis_,false);
    }
    Iterator end() const
    {
      return Iterator(mBasis_,true);
    }
    // for debuging
    int size() const
    {
      return mBasis_.size()*d;
    }
    operator F*()
    {
      return &(mBasis_[0]);
    }
    std::vector<F> mBasis_;
  };

  template< int dim, int dimR, class SF, class CF >
  struct VecLagrangeBasisCreator
  {
    typedef VirtualMonomialBasis<dim,SF> MBasis;
    typedef SF StorageField;
    typedef AlgLib::MultiPrecision< Precision<CF>::value > ComputeField;
    static const int dimension = dim;
    typedef PolynomialBasisWithMatrix<MBasis,StorageField,dimR,VecContainer<dimR,StorageField> > Basis;
    typedef unsigned int Key;

    template <class Topology>
    struct Maker
    {
      static void apply(unsigned int order,Basis* &basis)
      {
        const MBasis &virtBasis = MonomialBasisProvider<dimension,StorageField>::template basis<Topology>(order);
        basis = new Basis(virtBasis,order);
        VecLagrangeMatrix<Topology,ComputeField,dimR> matrix(order);
        basis->fill(matrix);
        std::stringstream name;
        name << "lagrange_" << Topology::name() << "_p" << order;
        basis->template printBasis<Topology>(name.str(),matrix);
      }
    };
  };

  template< int dim, int dimR, class SF, class CF = typename ComputeField< SF, 512 >::Type >
  struct VecLagrangeBasisProvider
    : public BasisProvider<VecLagrangeBasisCreator<dim,dimR,SF,CF> >
  {};
}
#endif // DUNE_ORTHONORMALBASIS_HH
