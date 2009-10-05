// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_BASISEVALUATOR_HH
#define DUNE_BASISEVALUATOR_HH

#include <vector>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <dune/grid/genericgeometry/topologytypes.hh>
#include <dune/grid/genericgeometry/misc.hh>

#include <dune/common/field.hh>

#include <dune/finiteelements/multiindex.hh>

namespace Dune
{
  template <class F,int dimD,int dimR,unsigned int deriv,bool useAll>
  struct Tensor
  {
    typedef Tensor<F,dimD,dimR,deriv-1,useAll> Base;
    static const int single = Base::single*dimD;
    static const int all = Base::all+single;
    typedef Dune::FieldVector<Dune::FieldVector<F,single>,dimR> Single;
    typedef Dune::FieldVector<Dune::FieldVector<F,all>,dimR> All;

    static const int blockSize = ( useAll ? all*dimR : single*dimR );
    typedef typename SelectType<useAll,All,Single>::Type Block;
    typedef Block Range;
  };
  template <class F,int dimD,int dimR,bool useAll>
  struct Tensor<F,dimD,dimR,2,useAll>
  {
    struct HJV
    {
      typedef FieldVector<FieldMatrix<F,dimD,dimD>,dimR> Hessian;
      typedef FieldMatrix<F,dimR,dimD> Jacobian;
      typedef FieldVector<F,dimR> Value;
      Hessian hessian;
      Jacobian jacobian;
      Value value;
    };
    typedef Tensor<F,dimD,dimR,1,useAll> Base;
    static const int single = Base::single*dimD;
    static const int all = Base::all+single;
    typedef Dune::FieldVector<Dune::FieldVector<F,single>,dimR> Single;
    typedef Dune::FieldVector<Dune::FieldVector<F,all>,dimR> All;

    static const int blockSize = ( useAll ? all*dimR : single*dimR );
    typedef typename SelectType<useAll,All,Single>::Type Block;
    typedef typename SelectType<useAll,HJV,typename HJV::Hessian>::Type Range;
  };
  template <class F,int dimD,int dimR,bool useAll>
  struct Tensor<F,dimD,dimR,1,useAll>
  {
    struct JV
    {
      typedef FieldMatrix<F,dimR,dimD> Jacobian;
      typedef FieldVector<F,dimR> Value;
      Jacobian jacobian;
      Value value;
    };
    typedef Tensor<F,dimD,dimR,0,useAll> Base;
    static const int single = Base::single*dimD;
    static const int all = Base::all+single;
    typedef Dune::FieldVector<Dune::FieldVector<F,single>,dimR> Single;
    typedef Dune::FieldVector<Dune::FieldVector<F,all>,dimR> All;

    static const int blockSize = ( useAll ? all*dimR : single*dimR );
    typedef typename SelectType<useAll,All,Single>::Type Block;
    typedef typename SelectType<useAll,JV,typename JV::Jacobian>::Type Range;
  };
  template <class F,int dimD,int dimR,bool useAll>
  struct Tensor<F,dimD,dimR,0,useAll>
  {
    struct V
    {
      typedef FieldVector<F,dimR> Value;
      Value value;
    };
    static const int single = 1;
    static const int all = 1;

    static const int blockSize = ( useAll ? all*dimR : single*dimR );
    typedef Dune::FieldVector<F,dimR> Block;
    typedef V Range;
  };

  template <class B>
  struct MonomialEvaluator
  {
    typedef B Basis;
    typedef typename Basis::Field Field;
    typedef typename Basis::DomainVector DomainVector;
    typedef std::vector<Field> Container;
    static const int dimension = Basis::dimension;
    template <class Tensor>
    struct BaseIterator
    {
      typedef typename Tensor::Block Block;
      typedef typename Tensor::Range Range;
      static const int blockSize = Tensor::blockSize;
      typedef typename Container::iterator CIter;
      BaseIterator(Container &container)
        : pos_(container.begin()), end_(container.end())
      {}
      const Range &operator*() const
      {
        assert(!done());
        return reinterpret_cast<const Range&>(*pos_);
      }
      Range &operator*()
      {
        assert(!done());
        return reinterpret_cast<Range&>(*pos_);
      }
      const Range *operator->() const
      {
        return &(operator*());
      }
      Range *operator->()
      {
        return &(operator*());
      }
      const Block block()
      {
        assert(!done());
        return reinterpret_cast<const Block&>(*pos_);
      }
      Block block() const
      {
        assert(!done());
        return reinterpret_cast<Block&>(*pos_);
      }
      bool done() const
      {
        return pos_==end_;
      }
      BaseIterator &operator++()
      {
        pos_ += blockSize;
        return *this;
      }
      BaseIterator &operator+=(unsigned int skip)
      {
        pos_ += skip*blockSize;
        return *this;
      }
    private:
      CIter pos_;
      const CIter end_;
    };
    /*
       MonomialEvaluator(const Basis &basis,unsigned int order)
       : basis_(basis),
       order_(order),
       size_(basis.size(order)),
       container_(0)
       {
       resize<2,true>();
       }
     */
    template <unsigned int deriv>
    struct Iterator
    {
      typedef BaseIterator<Tensor<Field,dimension,1,deriv,true> > All;
      typedef BaseIterator<Tensor<Field,dimension,1,deriv,false> > Single;
    };
    template <unsigned int deriv>
    typename Iterator<deriv>::Single evaluate(const DomainVector &x)
    {
      resize<deriv,false>();
      basis_.evaluate(order_,x,container_);
      return typename Iterator<deriv>::Single(container_);
    }
    template <unsigned int deriv>
    typename Iterator<deriv>::All evaluateAll(const DomainVector &x)
    {
      resize<deriv,true>();
      basis_.evaluate(order_,x,container_);
      return typename Iterator<deriv>::All(container_);
    }

    typename Iterator<0>::Single evaluate(const DomainVector &x)
    {
      return evaluate<0>(x);
    }
    typename Iterator<1>::Single jacobian(const DomainVector &x)
    {
      return evaluate<1>(x);
    }
    typename Iterator<0>::All evaluateAll(const DomainVector &x)
    {
      return evaluateAll<1>(x);
    }
    typename Iterator<1>::All jacobianAll(const DomainVector &x)
    {
      return evaluateAll<1>(x);
    }
    unsigned int order() const
    {
      return order_;
    }
    unsigned int size() const
    {
      return size_;
    }
  protected:
    MonomialEvaluator(const Basis &basis,unsigned int order,unsigned int size)
      : basis_(basis),
        order_(order),
        size_(size),
        container_(0)
    {
      resize<2,true>();
    }
    template <int deriv,bool useAll>
    void resize()
    {
      const int totalSize = Tensor<Field,dimension,1,deriv,useAll>::blockSize*size_;
      container_.resize(totalSize);
    }
    MonomialEvaluator(const MonomialEvaluator&);
    const Basis &basis_;
    unsigned int order_,size_;
    Container container_;
  };

  template <class B>
  struct StandardEvaluator : public MonomialEvaluator<B>
  {
    typedef B Basis;
    typedef typename Basis::Field Field;
    typedef typename Basis::DomainVector DomainVector;
    typedef std::vector<Field> Container;
    static const int dimension = Basis::dimension;
    typedef MonomialEvaluator<B> Base;

    template <unsigned int deriv>
    struct Iterator : public Base::template Iterator<deriv>
    {};

    StandardEvaluator(const Basis &basis,unsigned int order)
      : Base(basis,order,basis.size())
    {}
    template <unsigned int deriv>
    typename Iterator<deriv>::Single evaluate(const DomainVector &x)
    {
      Base::template resize<deriv,false>();
      basis_.evaluate(x,container_);
      return typename Iterator<deriv>::Single(container_);
    }
    template <unsigned int deriv>
    typename Iterator<deriv>::All evaluateAll(const DomainVector &x)
    {
      this->template resize<deriv,true>();
      basis_.evaluate(x,container_);
      return typename Iterator<deriv>::All(container_);
    }

    typename Iterator<0>::Single evaluate(const DomainVector &x)
    {
      return evaluate<0>(x);
    }
    typename Iterator<1>::Single jacobian(const DomainVector &x)
    {
      return evaluate<0>(x);
    }
    typename Iterator<0>::All evaluateAll(const DomainVector &x)
    {
      return evaluate<0>(x);
    }
    typename Iterator<1>::All jacobianAll(const DomainVector &x)
    {
      return evaluate<0>(x);
    }
  private:
    StandardEvaluator(const StandardEvaluator&);
    using Base::basis_;
    using Base::container_;
  };

  template <class B,class Fill>
  struct VecEvaluator : public StandardEvaluator<B>
  {
    typedef B Basis;
    typedef typename Basis::Field Field;
    typedef typename Basis::DomainVector DomainVector;
    typedef std::vector<Field> Container;
    static const int dimension = Basis::dimension;
    typedef StandardEvaluator<B> Base;
    static const int dimRange = Fill::dimRange;

    template <unsigned int deriv>
    struct Iterator
    {
      typedef typename Base::template BaseIterator<Tensor<Field,dimension,dimRange,deriv,true> > All;
      typedef typename Base::template BaseIterator<Tensor<Field,dimension,dimRange,deriv,false> > Single;
    };

    VecEvaluator(const Basis &basis,
                 unsigned int order,
                 const Fill &fill)
      : Base(basis,order),
        fill_(fill)
    {
      resize<2,true>();
    }
    template <unsigned int deriv>
    typename Iterator<deriv>::Single evaluate(const DomainVector &x)
    {
      fill_( x,Base::template evaluate<deriv>(x), vecContainer_ );
      return typename Iterator<deriv>::Single(vecContainer_);
    }
    template <unsigned int deriv>
    typename Iterator<deriv>::All evaluateAll(const DomainVector &x)
    {
      fill_( x,Base::template evaluateAll<deriv>(x), vecContainer_ );
      return typename Iterator<deriv>::All(vecContainer_);
    }
    typename Iterator<0>::Single evaluate(const DomainVector &x)
    {
      return evaluate<0>(x);
    }
    typename Iterator<1>::Single jacobian(const DomainVector &x)
    {
      return evaluate<1>(x);
    }
    typename Iterator<0>::All evaluateAll(const DomainVector &x)
    {
      return evaluateAll<0>(x);
    }
    typename Iterator<1>::All jacobianAll(const DomainVector &x)
    {
      return evaluateAll<1>(x);
    }
    unsigned int size() const
    {
      return size_*dimRange;
    }
  protected:
    template <int deriv,bool useAll>
    void resize()
    {
      const int totalSize = Tensor<Field,dimension,dimRange,deriv,useAll>::blockSize*size_;
      vecContainer_.resize(totalSize);
    }
    VecEvaluator(const VecEvaluator&);
    Container vecContainer_;
    using Base::size_;
    const Fill &fill_;
  };

  template <int dimR>
  struct DiagonalFill
  {
    static const int dimRange = dimR;
    template <class Domain, class Iter,class Field>
    void operator()(const Domain &x,
                    Iter iter,std::vector<Field> &vecContainer) const
    {
      typedef std::vector<Field> Container;
      typename Container::iterator vecIter = vecContainer.begin();
      for ( ; !iter.done(); ++iter)
      {
        const typename Iter::Block &block = iter.block();
        for (int r1=0; r1<dimR; ++r1)
        {
          for (int b=0; b<iter.blockSize; ++b)
          {
            for (int r2=0; r2<dimR; ++r2)
            {
              *vecIter = (r1==r2 ? block[b] : Field(0));
              ++vecIter;
            }
          }
        }
      }
    }
  };

  template <class B,int dimR>
  struct VectorialEvaluator
    : public VecEvaluator<B,DiagonalFill<dimR> >
  {
    typedef DiagonalFill<dimR> Fill;
    typedef VecEvaluator< B,Fill > Base;
    VectorialEvaluator(const B &basis,
                       unsigned int order)
      : Base(basis,order,fill_)
    {}
  private:
    Fill fill_;
  };

#if 0
  template <class B,class F>
  struct MultiIndexEvaluator
  {
    typedef B Basis;
    typedef F Field;
    typedef std::vector<typename Basis::Field> Container;
    static const int dimension = Basis::dimension;
    typedef Dune::FieldVector<Field,dimension> DomainVector;
    template <class BlockType,int deriv>
    struct BaseIterator
    {
      static const int blockSize = sizeof(BlockType)/sizeof(Field);
      typedef BlockType RangeVector;
      typedef typename Container::const_iterator CIter;
      BaseIterator(const std::vector<DomainVector> &x,
                   const Container &container)
        : pos_(container.begin()), end_(container.end()),
          x_(x)
      {
        set();
      }
      const RangeVector &operator*() const
      {
        assert(!done());
        return reinterpret_cast<const RangeVector&>(val_[0]);
      }
      const RangeVector *operator->() const
      {
        return &(operator*());
      }
      bool done() const
      {
        return pos_==end_;
      }
      BaseIterator operator++()
      {
        pos_ += 1;
        if (!done())
          set();
        return *this;
      }
      BaseIterator &operator+=(unsigned int skip)
      {
        pos_ += skip;
        if (!done())
          set();
        return *this;
      }
    private:
      void set()
      {
        if (deriv==0)
        {
          val_[0]=1;
          for (int d=0; d<dimension; ++d)
          {
            unsigned int o = pos_->z(d);
            assert( o<x_.size() );
            val_[0]  *= x_[o][d];
          }
        }
        else if (deriv==1)
        {
          for (int i=0; i<dimension; ++i)
          {
            unsigned int o = pos_->z(i);
            if ( o == 0)
              val_[i] = 0.;
            else
            {
              val_[i] = o;
              for (int d=0; d<dimension; ++d)
              {
                unsigned int o = pos_->z(d);
                o -= (d==i);
                assert( o<x_.size() );
                val_[i]  *= x_[o][d];
              }
            }
          }
          if (blockSize>dimension || deriv==0)
          {
            val_[dimension]=1;
            for (int d=0; d<dimension; ++d)
            {
              unsigned int o = pos_->z(d);
              assert( o<x_.size() );
              val_[dimension]  *= x_[o][d];
            }
          }
        }
      }
      CIter pos_;
      const CIter end_;
      const std::vector<DomainVector> &x_;
      Field val_[blockSize];
    };
    MultiIndexEvaluator(const Basis &basis,unsigned int order)
      : basis_(basis),
        order_(order),
        container_(basis.size(order)),
        x_(order+1)
    {
      typename Basis::DomainVector x;
      for( int i = 0; i < dimension; ++i )
        x[ i ].set( i, 1 );
      basis_.evaluate(order_,x,container_);
    }
    template <unsigned int deriv>
    struct Iterator
    {
      typedef BaseIterator<typename Tensor<Field,dimension,1,deriv,all>,deriv> All;
      typedef BaseIterator<typename Tensor<Field,dimension,1,deriv,all>,deriv> Single;
    };
    typename Iterator<0>::Single evaluate(const DomainVector &x)
    {
      setX(x);
      return typename Iterator<0>::Single(x_,container_);
    }
    template <unsigned int deriv>
    typename Iterator<deriv>::Single evaluate(const DomainVector &x)
    {
      setX(x);
      return typename Iterator<deriv>::Single(x_,container_);
    }
    template <unsigned int deriv>
    typename Iterator<deriv>::All evaluateAll(const DomainVector &x)
    {
      setX(x);
      return typename Iterator<deriv>::All(x_,container_);
    }
    typename Iterator<1>::Single jacobian(const DomainVector &x)
    {
      setX(x);
      return typename Iterator<1>::Single(x_,container_);
    }
    void setX(const DomainVector &x)
    {
      for (int d=0; d<dimension; ++d)
      {
        x_[0][d] = 1;
        for (unsigned int i=1; i<=order_; ++i) {
          x_[i][d]=x_[i-1][d]*x[d];
        }
      }
    }
    unsigned int order() const
    {
      return order_;
    }
    unsigned int size() const
    {
      return basis_.size(order);
    }
  private:
    MultiIndexEvaluator(const MultiIndexEvaluator&);
    const Basis &basis_;
    unsigned int order_;
    Container container_;
    std::vector<DomainVector> x_;
  };
#endif
}

#endif
