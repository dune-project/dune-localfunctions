// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
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
      typedef BaseIterator<Derivatives<Field,dimension,dimRange,deriv,DerivativeLayoutNS::derivative> > All;
      typedef BaseIterator<Derivatives<Field,dimension,1,0,DerivativeLayoutNS::value> > Integrate;
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
      const int totalSize = Derivatives<Field,dimension,dimRange,deriv,DerivativeLayoutNS::derivative>::size*size_;
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
    static const DerivativeLayoutNS::DerivativeLayout layout = Deriv::layout;
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

}

#endif
