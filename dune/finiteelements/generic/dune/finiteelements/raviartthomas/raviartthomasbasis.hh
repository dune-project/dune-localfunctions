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
  template <int d,class F>
  struct RaviartThomasContainer
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
    RaviartThomasContainer(int size)
      : mBasis_(size/d) {}
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

  template <class Topology,class scalar_t>
  struct RaviartThomasMatrix {
    enum {dim = Topology::dimension};
    typedef Dune::AlgLib::Matrix< scalar_t > mat_t;
    typedef StandardMonomialBasis<dim,scalar_t> MBasis;
    RaviartThomasMatrix(int order) : basis_(), order_(order)
    {}
    int colSize(int row) const {
      return (dim+1)*basis_.size(order_)-basis_.sizes(order_)[order_-1];
    }
    int rowSize() const {
      return (dim+1)*basis_.size(order_)-basis_.sizes(order_)[order_-1];
    }
    const Dune::FieldMatrix< scalar_t,dim,dim > operator() ( int r, int c ) const
    {
      Dune::FieldMatrix< scalar_t, dim,dim > ret(0);
      for (int i=0; i<dim; ++i)
        ret[i][i] = 1.;
      return ret;
    }
    void print(std::ostream& out,int N = rowSize()) const {}
    MBasis basis_;
    int order_;
    mat_t matrix_;
  };

  template< int dim, class SF, class CF >
  struct RaviartThomasBasisCreator
  {
    typedef StandardMonomialBasis<dim,SF> MBasis;
    typedef SF StorageField;
    typedef AlgLib::MultiPrecision< Precision<CF>::value > ComputeField;
    static const int dimension = dim;
    typedef PolynomialBasisWithMatrix<MBasis,StorageField,dimension,RaviartThomasContainer<dimension,StorageField> > Basis;
    typedef unsigned int Key;

    template <class Topology>
    struct Maker
    {
      static void apply(unsigned int order,Basis* &basis)
      {
        // bool RTBasis_only_for_implemented_for_simplex = GenericGeometry::IsSimplex<Topology>::value ;
        // assert(RTBasis_only_for_implemented_for_simplex);
        static MBasis _mBasis;
        basis = new Basis(_mBasis,order);
        RaviartThomasMatrix<Topology,ComputeField> matrix(order);
        basis->fill(matrix);
        std::stringstream name;
        name << "lagrange_" << Topology::name() << "_p" << order;
        basis->template printBasis<Topology>(name.str(),matrix);
      }
    };
  };

  template< int dim, class SF, class CF = typename ComputeField< SF, 512 >::Type >
  struct RaviartThomasBasisProvider
    : public BasisProvider<RaviartThomasBasisCreator<dim,SF,CF> >
  {};
}
#endif // DUNE_ORTHONORMALBASIS_HH
