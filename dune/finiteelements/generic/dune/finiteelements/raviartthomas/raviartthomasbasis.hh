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
      Iterator(const std::vector<F>& mBasis)
        : mBasis_(mBasis), pos_(0), dim_(0) {
        val_[dim_] = mBasis_[pos_];
      }
      Iterator(const Iterator& other)
        : mBasis_(other.mBasis_), pos_(other.pos_), dim_(other.dim_)
      {}
      Iterator &operator=(const Iterator& other)
      {
        mBasis_ = other.mBasis;
        pos_ = other.pos_;
        dim_ = other.dim_;
      }
      const FieldVector<F,d>& operator*() const
      {
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
        val_[dim_] = mBasis_[pos_];
        return *this;
      }
      FieldVector<F,d> val_;
      const std::vector<F>& mBasis_;
      int pos_,dim_;
    };
    typedef Dune::FieldVector<F,d> value_type;
    typedef Iterator const_iterator;
    RaviartThomasContainer(int size)
      : mBasis_(size) {}
    Iterator begin() const
    {
      return Iterator(mBasis_);
    }
    std::vector<F> mBasis_;
  };
  template <int dim,class SF>
  struct RaviartThomasMonomialBasis
  {
    typedef VirtualMonomialBasis<dim,SF> MBasis;
    static const int dimension = dim;
    typedef typename MBasis :: DomainVector DomainVector;
    RaviartThomasMonomialBasis(const VirtualMonomialBasis<dim,SF>& basis)
      : basis_(basis) {}
    int size(int order) const
    {
      return basis_.size(order);
    }
    void evaluate(unsigned int order, const DomainVector &x,
                  RaviartThomasContainer<dimension,SF> &val ) const
    {
      basis_.evaluate(order,x,val.mBasis_);
    }
    const VirtualMonomialBasis<dimension,SF>& basis_;
  };
  template< int dim, class SF >
  struct RaviartThomasMonomialBasisCreator
  {
    typedef RaviartThomasMonomialBasis<dim,SF> Basis;
    typedef SF StorageField;
    static const int dimension = dim;
    typedef unsigned int Key;

    template <class Topology>
    struct Maker
    {
      static void apply(unsigned int order,Basis* &basis)
      {
        basis = new Basis(MonomialBasisProvider<dimension,StorageField>::template basis<Topology>(order));
      }
    };
  };
  template< int dim, class SF>
  struct RaviartThomasMonomialBasisProvider
    : public BasisProvider<RaviartThomasMonomialBasisCreator<dim,SF> >
  {};
  // &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

  template <class Topology,class scalar_t>
  struct RaviartThomasMatrix {
    enum {dim = Topology::dimension};
    typedef Dune::AlgLib::Matrix< scalar_t > mat_t;
    RaviartThomasMatrix(int order)
    {}
    int colSize(int row) const {
      return matrix_.cols()/dim;
    }
    int rowSize() const {
      return matrix_.rows()/dim;
    }
    const Dune::FieldMatrix< scalar_t,dim,dim > operator() ( int r, int c ) const
    {
      Dune::FieldMatrix< scalar_t, dim,dim > ret(0);
      for (int i=0; i<dim; ++i)
        for (int j=0; j<dim; ++j)
          ret[i][j] = matrix_(c*dim+i,r*dim+j);
      return ret;
    }
    void print(std::ostream& out,int N = rowSize()) const {}
    mat_t matrix_;
  };


  template< int dim, class SF, class CF >
  struct RaviartThomasBasisCreator
  {
    typedef RaviartThomasMonomialBasis<dim,SF> RTBasis;
    typedef SF StorageField;
    typedef AlgLib::MultiPrecision< Precision<CF>::value > ComputeField;
    static const int dimension = dim;
    typedef PolynomialBasisWithMatrix<RTBasis,StorageField,dimension,RaviartThomasContainer<dimension,StorageField> > Basis;
    typedef unsigned int Key;

    template <class Topology>
    struct Maker
    {
      static void apply(unsigned int order,Basis* &basis)
      {
        const RTBasis &vBasis = RaviartThomasMonomialBasisProvider<dimension,StorageField>::template basis<Topology>(order);
        basis = new Basis(vBasis,order);
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
