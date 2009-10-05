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
  template <class MBasis>
  struct RaviartThomasFill
  {
    static const int dimRange = MBasis::dimension;
    RaviartThomasFill(const MBasis &basis)
      : basis_(basis)
    {}
    template <class Domain,class Iter,class Field>
    void operator()(const Domain &x, Iter iter,std::vector<Field> &vecContainer) const
    {
      typedef std::vector<Field> Container;
      typename Container::iterator vecIter = vecContainer.begin();
      unsigned int notHomogen = basis_.size(basis_.order()-1);
      std::cout << " Order: " << basis_.order()
                << " notHomo: " << notHomogen
                << " Homo:    " << basis_.size()-notHomogen
                << std::endl;
      for (unsigned int baseFunc = 0 ;
           baseFunc<notHomogen; ++iter, ++baseFunc)
      {
        const typename Iter::Block &block = iter.block();
        for (int b=0; b<iter.blockSize; ++b)
        {
          for (int r1=0; r1<dimRange; ++r1)
          {
            for (int r2=0; r2<dimRange; ++r2)
            {
              *vecIter = (r1==r2 ? block[b] : Field(0));
              ++vecIter;
            }
          }
        }
      }
      for ( ; !iter.done(); ++iter )
      {
        const typename Iter::Block &block = iter.block();
        for (int b=0; b<iter.blockSize; ++b)
        {
          for (int r1=0; r1<dimRange; ++r1)
          {
            for (int r2=0; r2<dimRange; ++r2)
            {
              *vecIter = (r1==r2 ? block[b] : Field(0));
              ++vecIter;
            }
          }
        }
        for (int b=0; b<iter.blockSize; ++b)
        {
          for (int r1=0; r1<dimRange; ++r1)
          {
            *vecIter = x[r1]*block[b];
            ++vecIter;
          }
        }
      }
    }
    const MBasis &basis_;
  };

  template <class B>
  struct RaviartThomasEvaluator
    : public VecEvaluator<B,RaviartThomasFill<B> >
  {
    typedef RaviartThomasFill<B> Fill;
    typedef VecEvaluator< B,Fill > Base;
    RaviartThomasEvaluator(const B &basis,
                           unsigned int order)
      : Base(basis,order,fill_),
        fill_(basis)
    {}
  private:
    Fill fill_;
  };

  template <class Topology,class scalar_t>
  struct RaviartThomasMatrix {
    enum {dim = Topology::dimension};
    typedef Dune::AlgLib::Matrix< scalar_t > mat_t;
    typedef StandardMonomialBasis<dim,scalar_t> MBasis;
    RaviartThomasMatrix(int order) : basis_(order), order_(order)
    {}
    int colSize(int row) const {
      return (dim+1)*basis_.size(order_)-basis_.sizes(order_)[order_-1];
    }
    int rowSize() const {
      return (dim+1)*basis_.size(order_)-basis_.sizes(order_)[order_-1];
    }
    const scalar_t operator() ( int r, int c ) const
    {
      return ( (r==c) ? scalar_t(1) : scalar_t(0) );
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
    typedef PolynomialBasisWithMatrix<RaviartThomasEvaluator<MBasis>,SparseCoeffMatrix<StorageField> > Basis;
    typedef unsigned int Key;

    template <class Topology>
    struct Maker
    {
      static void apply(unsigned int order,Basis* &basis)
      {
        // bool RTBasis_only_for_implemented_for_simplex = GenericGeometry::IsSimplex<Topology>::value ;
        // assert(RTBasis_only_for_implemented_for_simplex);
        static MBasis _mBasis(order);
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
