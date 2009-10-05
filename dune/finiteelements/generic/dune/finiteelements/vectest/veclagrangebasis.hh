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
    const scalar_t operator() ( int r, int c ) const
    {
      if (c%dimR == r%dimR)
        return matrix_(c/dimR,r/dimR);
      else
        return Zero<scalar_t>();
    }
    void print(std::ostream& out,int N = rowSize()) const {}
    mat_t matrix_;
  };

  template <class Evaluator,int d>
  struct VecEvaluator
  {
    typedef Evaluator Container;
    typedef typename Evaluator::Basis Basis;
    typedef typename Evaluator::DomainVector DomainVector;
    typedef typename Evaluator::Field Field;
    static const int dimension = Evaluator :: dimension;

    struct Iterator
    {
      typedef typename Container::Iterator CIter;
      typedef Dune::FieldVector<typename CIter::RangeVector,d> RangeVector;
      Iterator(const CIter& pos)
        : pos_(pos),
          dim_(0)
      {
        val_[dim_] = *pos_;
      }
      /*
         Iterator(const Iterator& other)
         : val_(other.val_),
         eval_(other.eval_),
         pos_(other.pos_),
         dim_(other.dim_)
         {
         }
         Iterator &operator=(const Iterator& other)
         {
         val_ = other.val_;
         eval_ = other.eval_;
         pos_ = other.pos_;
         dim_ = other.dim_;
         }
       */
      bool done()
      {
        return pos_.done();
      }
      const RangeVector& operator*() const
      {
        assert(!pos_.done());
        return val_;
      }
      const Iterator &operator++()
      {
        val_[dim_] = 0;
        ++dim_;
        if (dim_ == d)
        {
          dim_ = 0;
          ++pos_;
          if (pos_.done())
            return *this;
        }
        val_[dim_] = *pos_;
        return *this;
      }
    private:
      typename Container::Iterator pos_;
      RangeVector val_;
      int dim_;
    };

    VecEvaluator(const Basis &basis, unsigned int order)
      : eval_(basis,order)
    {}
    Iterator evaluate(const DomainVector &x)
    {
      return Iterator(eval_.evaluate(x));
    }
    /*
       void evaluate(const DomainVector &x)
       {
       eval_.evaluate(x);
       }
       Iterator begin() const
       {
       return Iterator(eval_,false);
       }
       Iterator end() const
       {
       return Iterator(eval_,true);
       }
     */
    unsigned int order() const
    {
      return eval_.order();
    }
    unsigned int size() const
    {
      return eval_.size()*d;
    }
  private:
    Evaluator eval_;
  };

  template< int dim, int dimR, class SF, class CF >
  struct VecLagrangeBasisCreator
  {
    typedef VirtualMonomialBasis<dim,SF> MBasis;
    typedef SF StorageField;
    typedef AlgLib::MultiPrecision< Precision<CF>::value > ComputeField;
    static const int dimension = dim;
    typedef VecEvaluator<StandardEvaluator<MBasis>,dimR> Evaluator;
    typedef PolynomialBasisWithMatrix<Evaluator> Basis;
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
