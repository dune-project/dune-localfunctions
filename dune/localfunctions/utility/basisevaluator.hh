// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_BASISEVALUATOR_HH
#define DUNE_BASISEVALUATOR_HH

#include <vector>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/typetraits.hh>

#include <dune/geometry/type.hh>

#include <dune/localfunctions/utility/field.hh>
#include <dune/localfunctions/utility/multiindex.hh>
#include <dune/localfunctions/utility/tensor.hh>

namespace Dune
{
  /*******************************************
  * Should be removed as soon as the Tensor
  * classes have been revisited. See remarks
  * in tensor.hh (also hold true here).
  *******************************************/


  template <class B>
  struct MonomialEvaluator
  {
    typedef B Basis;
    typedef typename Basis::Field Field;
    typedef typename Basis::DomainVector DomainVector;
    static const int dimension = Basis::dimension;
    static const int dimRange = Basis::dimRange;

    typedef std::vector<Field> Container;

    template< class Deriv >
    struct BaseIterator;

    template <unsigned int deriv>
    struct Iterator
    {
      typedef BaseIterator<Derivatives<Field,dimension,dimRange,deriv,derivative> > All;
      typedef BaseIterator<Derivatives<Field,dimension,1,0,value> > Integrate;
    };

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
    {}
    template <int deriv>
    void resize()
    {
      const int totalSize = Derivatives<Field,dimension,dimRange,deriv,derivative>::size*size_;
      container_.resize(totalSize);
    }
    MonomialEvaluator(const MonomialEvaluator&);
    const Basis &basis_;
    unsigned int order_,size_;
    Container container_;
  };


  template< class B >
  template< class Deriv >
  struct MonomialEvaluator< B >::BaseIterator
  {
    typedef Deriv Derivatives;
    typedef typename Deriv::Field Field;
    static const unsigned int blockSize = Deriv::size;
    typedef Dune::FieldVector<Field,blockSize> Block;
    static const DerivativeLayout layout = Deriv::layout;
    static const unsigned int dimDomain = Deriv::dimDomain;
    static const unsigned int dimRange = Deriv::dimRange;

    typedef std::vector<Field> Container;
    typedef typename Container::iterator CIter;

    explicit BaseIterator ( Container &container )
      : pos_( container.begin() ),
        end_( container.end() )
    {}

    const Deriv &operator*() const
    {
      assert(!done());
      return reinterpret_cast<const Deriv&>(*pos_);
    }

    const Deriv *operator->() const
    {
      return &(operator*());
    }

    bool done () const
    {
      return pos_ == end_;
    }

    BaseIterator &operator++ ()
    {
      pos_ += blockSize;
      return *this;
    }

    BaseIterator &operator+= ( unsigned int skip )
    {
      pos_ += skip*blockSize;
      return *this;
    }

  private:
    CIter pos_;
    const CIter end_;
  };

  template< class B >
  struct StandardEvaluator
    : public MonomialEvaluator< B >
  {
    typedef B Basis;
    typedef typename Basis::Field Field;
    typedef typename Basis::DomainVector DomainVector;
    typedef std::vector<Field> Container;
    static const int dimension = Basis::dimension;
    static const int dimRange = Basis::dimRange;
    typedef MonomialEvaluator<B> Base;

    template <unsigned int deriv>
    struct Iterator : public Base::template Iterator<deriv>
    {};

    StandardEvaluator(const Basis &basis)
      : Base(basis,basis.order(),basis.size())
    {}
    template <unsigned int deriv,class DVector>
    typename Iterator<deriv>::All evaluate(const DVector &x)
    {
      Base::template resize<deriv>();
      basis_.template evaluate<deriv>(x,&(container_[0]));
      return typename Iterator<deriv>::All(container_);
    }
    typename Iterator<0>::Integrate integrate()
    {
      Base::template resize<0>();
      basis_.integrate(&(container_[0]));
      return typename Iterator<0>::Integrate(container_);
    }

  protected:
    StandardEvaluator ( const Basis &basis, unsigned int size )
      : Base( basis, basis.order(), size )
    {}

  private:
    StandardEvaluator(const StandardEvaluator&);
    using Base::basis_;
    using Base::container_;
  };

#if 0 // OLD OLD
  template< class B, class Fill >
  struct VecEvaluator
    : public StandardEvaluator< B >
  {
    typedef B Basis;
    typedef typename Basis::Field Field;
    static const int dimension = Basis::dimension;
    static const int dimRange = Basis::dimRange*Fill::dimRange;
    typedef typename Basis::DomainVector DomainVector;
    typedef std::vector<Field> Container;
    typedef StandardEvaluator<B> Base;

    template <unsigned int deriv>
    struct Iterator
    {
      typedef typename Base::template BaseIterator<Derivatives<Field,dimension,dimRange,deriv,Fill::layout> > All;
    };

    VecEvaluator ( const Basis &basis, const Fill &fill )
      : Base( basis, basis.size() ),
        fill_( fill ),
        size_( basis.size()*dimRange )
    {}
    template <unsigned int deriv>
    typename Iterator<deriv>::All evaluate(const DomainVector &x)
    {
      resize< deriv >();
      fill_.template apply<deriv>( x,Base::template evaluate<deriv>(x), vecContainer_ );
      std::vector<Derivatives<Field,dimension,dimRange,deriv,Fill::layout> >& derivContainer =
        reinterpret_cast<std::vector<Derivatives<Field,dimension,dimRange,deriv,Fill::layout> >&>(vecContainer_);
      return typename Iterator<deriv>::All(derivContainer);
    }
    template <unsigned int deriv,class DVector>
    typename Iterator<deriv>::All evaluate(const DVector &x)
    {
      resize< deriv >();
      fill_.template apply<deriv>( x,Base::template evaluate<deriv>(x), vecContainer_ );
      std::vector<Derivatives<Field,dimension,dimRange,deriv,Fill::layout> >& derivContainer =
        reinterpret_cast<std::vector<Derivatives<Field,dimension,dimRange,deriv,Fill::layout> >&>(vecContainer_);
      return typename Iterator<deriv>::All(derivContainer);
    }
    unsigned int size() const
    {
      return size_;
    }

  protected:
    VecEvaluator ( const Basis &basis, const Fill &fill, unsigned int size )
      : Base( basis, basis.size() ),
        fill_( fill ),
        size_( size )
    {
      resize< 2 >();
    }

    template <int deriv>
    void resize()
    {
      const int totalSize = Derivatives<Field,dimension,dimRange,deriv,derivative>::size*size_;
      vecContainer_.resize(totalSize);
    }

    VecEvaluator(const VecEvaluator&);

    Container vecContainer_;
    const Fill &fill_;
    unsigned int size_;
  };

  template <int dimR,DerivativeLayout layout>
  struct DiagonalFill;

  template <int dimR>
  struct DiagonalFill<dimR,derivative>
  {
    static const DerivativeLayout layout = derivative;
    static const int dimRange = dimR;
    template <int deriv, class Domain, class Iter,class Field>
    void apply(const Domain &x,
               Iter iter,std::vector<Field> &vecContainer) const
    {
      typedef std::vector<Field> Container;
      typename Container::iterator vecIter = vecContainer.begin();
      for ( ; !iter.done(); ++iter)
      {
        const typename Iter::Block &block = iter->block();
        for (int r1=0; r1<dimR; ++r1)
        {
          unsigned int b = 0;
          apply<Field>(r1,x,block,b,vecIter);
        }
      }
    }
    template <class Field, class Domain, class Block,class VecIter>
    void apply(int r1, const Domain &x,
               const Block &block,unsigned int &b,
               VecIter &vecIter) const
    {
      unsigned int bStart = b;
      unsigned int bEnd = b+Block::size;
      apply<Field>(r1,x,block,bStart,bEnd,vecIter);
      b=bEnd;
    }
    template <class Field, class Domain, class Block,class VecIter>
    void apply(int r1, const Domain &x,const Block &block,
               unsigned int bStart, unsigned int bEnd,
               VecIter &vecIter) const
    {
      for (int r2=0; r2<dimR; ++r2)
      {
        for (unsigned int bb=bStart; bb<bEnd; ++bb)
        {
          *vecIter = (r1==r2 ? block[bb] : Field(0));
          ++vecIter;
        }
      }
    }
  };
  template <int dimR>
  struct DiagonalFill<dimR,value>
  {
    static const DerivativeLayout layout = value;
    static const int dimRange = dimR;
    template <int deriv, class Domain, class Iter,class Field>
    void apply(const Domain &x,
               Iter iter,std::vector<Field> &vecContainer) const
    {
      typedef std::vector<Field> Container;
      typename Container::iterator vecIter = vecContainer.begin();
      for ( ; !iter.done(); ++iter)
      {
        const typename Iter::Block &block = iter->block();
        for (int r1=0; r1<dimR; ++r1)
        {
          unsigned int b = 0;
          apply<Field>(std::integral_constant<int,deriv>(),r1,x,block,b,vecIter);
        }
      }
    }
    template <class Field, class Domain, class Block,class VecIter,int deriv>
    void apply(const integral_constat<int,deriv>&, int r1, const Domain &x,
               const Block &block,unsigned int &b,
               VecIter &vecIter) const
    {
      apply<Field>(std::integral_constant<int,deriv-1>(),r1,x,block,b,vecIter);
      unsigned int bStart = b;
      unsigned int bEnd = b+LFETensor<Field,Domain::dimension,deriv>::size;
      apply<Field>(r1,x,block,bStart,bEnd,vecIter);
      b=bEnd;
    }
    template <class Field, class Domain, class Block,class VecIter>
    void apply(const std::integral_constant<int,0>&, int r1, const Domain &x,
               const Block &block,unsigned int &b,
               VecIter &vecIter) const
    {
      apply<Field>(r1,x,block,b,b+1,vecIter);
      ++b;
    }
    template <class Field, class Domain, class Block,class VecIter>
    void apply(int r1, const Domain &x,const Block &block,
               unsigned int bStart, unsigned int bEnd,
               VecIter &vecIter) const
    {
      for (int r2=0; r2<dimR; ++r2)
      {
        for (unsigned int bb=bStart; bb<bEnd; ++bb)
        {
          *vecIter = (r1==r2 ? block[bb] : Field(0));
          ++vecIter;
        }
      }
    }
  };

  template <class B,int dimR,DerivativeLayout layout>
  struct VectorialEvaluator
    : public VecEvaluator<B,DiagonalFill<dimR,layout> >
  {
    typedef DiagonalFill<dimR,layout> Fill;
    typedef VecEvaluator< B,Fill > Base;
    VectorialEvaluator(const B &basis)
      : Base(basis,fill_,basis.size()*dimR)
    {}
  private:
    Fill fill_;
  };
#endif // OLD OLD

}

#endif
