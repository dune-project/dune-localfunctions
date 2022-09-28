// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

#ifndef DUNE_TENSOR_HH
#define DUNE_TENSOR_HH

#include <ostream>
#include <vector>

#include <dune/common/fvector.hh>

#include <dune/localfunctions/utility/field.hh>

namespace Dune
{
  /***********************************************
   * The classes here are work in progress.
   * Basically they provide tensor structures for
   * higher order derivatives of vector valued function.
   * Two storage structures are provided
   * (either based on the components of the vector valued
   * functions or on the order of the derivative).
   * Conversions are supplied between the two storage
   * structures and simple operations, which make the
   * code difficult to use and requires rewritting...
   ***************************************************/

  // Structure for scalar tensor of order deriv
  template <class F,int dimD,unsigned int deriv>
  class LFETensor
  {
    typedef LFETensor<F,dimD,deriv> This;
    typedef LFETensor<F,dimD-1,deriv> BaseDim;
    typedef LFETensor<F,dimD,deriv-1> BaseDeriv;

  public:
    typedef F field_type;
    static const unsigned int size = BaseDim::size+BaseDeriv::size;
    typedef Dune::FieldVector<F,size> Block;

    template< class FF >
    This &operator= ( const FF &f )
    {
      block() = field_cast< F >( f );
      return *this;
    }

    This &operator= ( const Block &b )
    {
      block() = b;
      return *this;
    }

    This &operator*= ( const field_type &f )
    {
      block() *= f;
      return *this;
    }

    const field_type &operator[] ( const unsigned int i ) const
    {
      return block()[ i ];
    }

    field_type &operator[] ( const unsigned int i )
    {
      return block()[ i ];
    }

    Block &block()
    {
      return block_;
    }
    const Block &block() const
    {
      return block_;
    }
    void axpy(const F& a, const This &y)
    {
      block().axpy(a,y.block());
    }
    template <class Fy>
    void assign(const LFETensor<Fy,dimD,deriv> &y)
    {
      field_cast(y.block(),block());
    }
    Block block_;
  };

  // ******************************************
  template <class F,unsigned int deriv>
  struct LFETensor<F,0,deriv>
  {
    static const int size = 0;
  };

  template <class F>
  struct LFETensor<F,0,0>
  {
    static const int size = 1;
  };

  template <class F,int dimD>
  class LFETensor<F,dimD,0>
  {
    typedef LFETensor<F,dimD,0> This;

  public:
    typedef F field_type;
    static const int size = 1;
    typedef Dune::FieldVector<F,size> Block;

    template< class FF >
    This &operator= ( const FF &f )
    {
      block() = field_cast< F >( f );
      return *this;
    }

    This &operator= ( const Block &b )
    {
      block() = b;
      return *this;
    }

    This &operator*= ( const field_type &f )
    {
      block() *= f;
      return *this;
    }

    const F &operator[] ( const unsigned int i ) const
    {
      return block()[ i ];
    }

    F &operator[] ( const unsigned int i )
    {
      return block()[ i ];
    }

    void axpy(const F& a, const This &y)
    {
      block().axpy(a,y.block());
    }
    template <class Fy>
    void assign(const LFETensor<Fy,dimD,0> &y)
    {
      field_cast(y.block(),block());
    }

    Block &block()
    {
      return block_;
    }
    const Block &block() const
    {
      return block_;
    }
    Block block_;
  };
  // ***********************************************************
  // Structure for all derivatives up to order deriv
  // for vector valued function
  namespace DerivativeLayoutNS {
    enum DerivativeLayout {value,derivative};
  }
  template <class F,int dimD,int dimR,unsigned int deriv,
      DerivativeLayoutNS::DerivativeLayout layout>
  struct Derivatives;

  // Implemnetation for valued based layout
  template <class F,int dimD,int dimR,unsigned int deriv>
  struct Derivatives<F,dimD,dimR,deriv,DerivativeLayoutNS::value>
    : public Derivatives<F,dimD,dimR,deriv-1,DerivativeLayoutNS::value>
  {
    typedef Derivatives<F,dimD,dimR,deriv,DerivativeLayoutNS::value> This;
    typedef Derivatives<F,dimD,dimR,deriv-1,DerivativeLayoutNS::value> Base;
    typedef LFETensor<F,dimD,deriv> ThisLFETensor;

    typedef F Field;
    typedef F field_type;

    static const DerivativeLayoutNS::DerivativeLayout layout = DerivativeLayoutNS::value;
    static const unsigned int dimDomain = dimD;
    static const unsigned int dimRange = dimR;
    constexpr static int size = Base::size+ThisLFETensor::size*dimR;
    typedef Dune::FieldVector<F,size> Block;

    This &operator=(const F& f)
    {
      block() = f;
      return *this;
    }
    This &operator=(const Dune::FieldVector<ThisLFETensor,dimR> &t)
    {
      tensor_ = t;
      return *this;
    }
    template <unsigned int dorder>
    This &operator=(const Dune::FieldVector<LFETensor<F,dimD,dorder>,dimR> &t)
    {
      tensor<dorder>() = t;
      return *this;
    }
    This &operator=(const Block &t)
    {
      block() = t;
      return *this;
    }

    This &operator*= ( const field_type &f )
    {
      block() *= f;
      return *this;
    }

    void axpy(const F &a, const This &y)
    {
      block().axpy(a,y.block());
    }

    // assign with same layout (only different Field)
    template <class Fy>
    void assign(const Derivatives<Fy,dimD,dimR,deriv,DerivativeLayoutNS::value> &y)
    {
      field_cast(y.block(),block());
    }
    // assign with different layout (same dimRange)
    template <class Fy>
    void assign(const Derivatives<Fy,dimD,dimR,deriv,DerivativeLayoutNS::derivative> &y)
    {
      Base::assign(y);
      for (int rr=0; rr<dimR; ++rr)
        tensor_[rr] = y[rr].template tensor<deriv>()[0];
    }
    // assign with rth component of function
    template <class Fy,int dimRy>
    void assign(const Derivatives<Fy,dimD,dimRy,deriv,DerivativeLayoutNS::value> &y,unsigned int r)
    {
      assign<Fy,dimRy>(y.block(),r);
    }
    // assign with scalar functions to component r
    template <class Fy>
    void assign(unsigned int r,const Derivatives<Fy,dimD,1,deriv,DerivativeLayoutNS::value> &y)
    {
      assign(r,y.block());
    }
    template <class Fy>
    void assign(unsigned int r,const Derivatives<Fy,dimD,1,deriv,DerivativeLayoutNS::derivative> &y)
    {
      assign(r,y[0]);
    }

    Block &block()
    {
      return reinterpret_cast<Block&>(*this);
    }
    const Block &block() const
    {
      return reinterpret_cast<const Block&>(*this);
    }

    template <unsigned int dorder>
    const Dune::FieldVector<LFETensor<F,dimD,dorder>,dimR> &tensor() const
    {
      // use integral_constant<int,...> here to stay compatible with Int2Type
      const std::integral_constant<int,dorder> a = {};
      return tensor(a);
    }
    template <unsigned int dorder>
    Dune::FieldVector<LFETensor<F,dimD,dorder>,dimR> &tensor()
    {
      // use integral_constant<int,...> here to stay compatible with Int2Type
      return tensor(std::integral_constant<int,dorder>());
    }
    template <unsigned int dorder>
    const Dune::FieldVector<F,LFETensor<F,dimD,dorder>::size*dimR> &block() const
    {
      // use integral_constant<int,...> here to stay compatible with Int2Type
      const std::integral_constant<int,dorder> a = {};
      return reinterpret_cast<const Dune::FieldVector<F,LFETensor<F,dimD,dorder>::size*dimR>&>(tensor(a));
    }
    template <unsigned int dorder>
    Dune::FieldVector<F,LFETensor<F,dimD,dorder>::size*dimR> &block()
    {
      // use integral_constant<int,...> here to stay compatible with Int2Type
      const std::integral_constant<int,dorder> a = {};
      return reinterpret_cast<Dune::FieldVector<F,LFETensor<F,dimD,dorder>::size*dimR>&>(tensor(a));
    }
    ThisLFETensor &operator[](int r) {
      return tensor_[r];
    }
    const ThisLFETensor &operator[](int r) const {
      return tensor_[r];
    }
  protected:
    template <class Fy,int dimRy>
    void assign(const FieldVector<Fy,size*dimRy> &y,unsigned int r)
    {
      Base::template assign<Fy,dimRy>(reinterpret_cast<const FieldVector<Fy,Base::size*dimRy>&>(y),r);
      tensor_[0] = reinterpret_cast<const FieldVector<Fy,ThisLFETensor::size>&>(y[Base::size*dimRy+r*ThisLFETensor::size]);
    }
    template <class Fy>
    void assign(unsigned int r,const FieldVector<Fy,size/dimR> &y)
    {
      Base::assign(r,reinterpret_cast<const FieldVector<Fy,Base::size/dimR>&>(y));
      tensor_[r] = reinterpret_cast<const FieldVector<Fy,ThisLFETensor::size>&>(y[Base::size/dimR]);
    }
    // assign with different layout (same dimRange)
    template <class Fy,unsigned int dy>
    void assign(const Derivatives<Fy,dimD,dimR,dy,DerivativeLayoutNS::derivative> &y)
    {
      Base::assign(y);
      for (int rr=0; rr<dimR; ++rr)
        tensor_[rr] = y[rr].template tensor<deriv>()[0];
    }

    template <int dorder>
    const Dune::FieldVector<LFETensor<F,dimD,dorder>,dimR> &
    tensor(const std::integral_constant<int,dorder> &dorderVar) const
    {
      return Base::tensor(dorderVar);
    }
    const Dune::FieldVector<LFETensor<F,dimD,deriv>,dimR> &
    tensor(const std::integral_constant<int,deriv> &dorderVar) const
    {
      return tensor_;
    }
    template <int dorder>
    Dune::FieldVector<LFETensor<F,dimD,dorder>,dimR> &
    tensor(const std::integral_constant<int,dorder> &dorderVar)
    {
      return Base::tensor(dorderVar);
    }
    Dune::FieldVector<LFETensor<F,dimD,deriv>,dimR> &
    tensor(const std::integral_constant<int,deriv> &dorderVar)
    {
      return tensor_;
    }
    Dune::FieldVector<ThisLFETensor,dimR> tensor_;
  };

  template <class F,int dimD,int dimR>
  struct Derivatives<F,dimD,dimR,0,DerivativeLayoutNS::value>
  {
    typedef Derivatives<F,dimD,dimR,0,DerivativeLayoutNS::value> This;
    typedef LFETensor<F,dimD,0> ThisLFETensor;

    typedef F Field;
    typedef F field_type;

    static const DerivativeLayoutNS::DerivativeLayout layout = DerivativeLayoutNS::value;
    static const unsigned int dimDomain = dimD;
    static const unsigned int dimRange = dimR;
    constexpr static int size = ThisLFETensor::size*dimR;
    typedef Dune::FieldVector<F,size> Block;

    template <class FF>
    This &operator=(const FF& f)
    {
      for (int r=0; r<dimR; ++r)
        tensor_[r] = field_cast<F>(f);
      return *this;
    }
    This &operator=(const Dune::FieldVector<ThisLFETensor,dimR> &t)
    {
      tensor_ = t;
      return *this;
    }

    This &operator=(const Block &t)
    {
      block() = t;
      return *this;
    }

    This &operator*= ( const field_type &f )
    {
      block() *= f;
      return *this;
    }

    void axpy(const F &a, const This &y)
    {
      block().axpy(a,y.block());
    }
    template <class Fy>
    void assign(const Derivatives<Fy,dimD,dimR,0,DerivativeLayoutNS::value> &y)
    {
      field_cast(y.block(),block());
    }
    template <class Fy>
    void assign(const Derivatives<Fy,dimD,dimR,0,DerivativeLayoutNS::derivative> &y)
    {
      for (int rr=0; rr<dimR; ++rr)
        tensor_[rr] = y[rr].template tensor<0>()[0];
    }
    template <class Fy,int dimRy>
    void assign(const Derivatives<Fy,dimD,dimRy,0,DerivativeLayoutNS::value> &y,unsigned int r)
    {
      assign<Fy,dimRy>(y.block(),r);
    }
    template <class Fy>
    void assign(unsigned int r,const Derivatives<Fy,dimD,1,0,DerivativeLayoutNS::value> &y)
    {
      tensor_[r].assign(y[0]);
    }
    template <class Fy>
    void assign(unsigned int r,const Derivatives<Fy,dimD,1,0,DerivativeLayoutNS::derivative> &y)
    {
      tensor_[r].assign(y[0][0]);
    }

    Block &block()
    {
      return reinterpret_cast<Block&>(*this);
    }
    const Block &block() const
    {
      return reinterpret_cast<const Block&>(*this);
    }

    ThisLFETensor &operator[](int r) {
      return tensor_[r];
    }
    const ThisLFETensor &operator[](int r) const {
      return tensor_[r];
    }
    template <int dorder>
    const Dune::FieldVector<LFETensor<F,dimD,0>,dimR> &tensor() const
    {
      return tensor_;
    }
    Dune::FieldVector<LFETensor<F,dimD,0>,dimR> &tensor()
    {
      return tensor_;
    }
    template <unsigned int dorder>
    const Dune::FieldVector<F,LFETensor<F,dimD,dorder>::size*dimR> &block() const
    {
      // use integral_constant<int,...> here to stay compatible with Int2Type
      const std::integral_constant<int,dorder> a = {};
      return reinterpret_cast<const Dune::FieldVector<F,LFETensor<F,dimD,dorder>::size*dimR>&>(tensor(a));
    }
    template <unsigned int dorder>
    Dune::FieldVector<F,LFETensor<F,dimD,dorder>::size*dimR> &block()
    {
      // use integral_constant<int,...> here to stay compatible with Int2Type
      const std::integral_constant<int,dorder> a = {};
      return reinterpret_cast<Dune::FieldVector<F,LFETensor<F,dimD,dorder>::size*dimR>&>(tensor(a));
    }

  protected:
    const Dune::FieldVector<LFETensor<F,dimD,0>,dimR> &
    tensor(const std::integral_constant<int,0> &dorderVar) const
    {
      return tensor_;
    }
    Dune::FieldVector<LFETensor<F,dimD,0>,dimR> &
    tensor(const std::integral_constant<int,0> &dorderVar)
    {
      return tensor_;
    }
    template <class Fy,unsigned int dy>
    void assign(const Derivatives<Fy,dimD,dimR,dy,DerivativeLayoutNS::derivative> &y)
    {
      for (int rr=0; rr<dimR; ++rr)
        tensor_[rr] = y[rr].template tensor<0>()[0];
    }
    template <class Fy,int dimRy>
    void assign(const FieldVector<Fy,size*dimRy> &y,unsigned int r)
    {
      tensor_[0] = reinterpret_cast<const FieldVector<Fy,ThisLFETensor::size>&>(y[r*ThisLFETensor::size]);
    }
    template <class Fy>
    void assign(unsigned int r,const FieldVector<Fy,size/dimR> &y)
    {
      tensor_[r] = y;
    }
    Dune::FieldVector<ThisLFETensor,dimR> tensor_;
  };

  // Implemnetation for DerivativeLayoutNS::derivative based layout
  template <class F,int dimD,int dimR,unsigned int deriv>
  struct Derivatives<F,dimD,dimR,deriv,DerivativeLayoutNS::derivative>
  {
    typedef Derivatives<F,dimD,dimR,deriv,DerivativeLayoutNS::derivative> This;
    typedef Derivatives<F,dimD,1,deriv,DerivativeLayoutNS::value> ScalarDeriv;

    typedef F Field;
    typedef F field_type;

    static const DerivativeLayoutNS::DerivativeLayout layout = DerivativeLayoutNS::value;
    static const unsigned int dimDomain = dimD;
    static const unsigned int dimRange = dimR;
    constexpr static int size = ScalarDeriv::size*dimR;
    typedef Dune::FieldVector<F,size> Block;

    template <class FF>
    This &operator=(const FF& f)
    {
      block() = field_cast<F>(f);
      return *this;
    }
    This &operator=(const Block &t)
    {
      block() = t;
      return *this;
    }

    This &operator*= ( const field_type &f )
    {
      block() *= f;
      return *this;
    }

    template <class FF>
    void axpy(const FF &a, const This &y)
    {
      block().axpy(field_cast<F>(a),y.block());
    }
    // assign with same layout (only different Field)
    template <class Fy>
    void assign(const Derivatives<Fy,dimD,dimR,deriv,DerivativeLayoutNS::derivative> &y)
    {
      field_cast(y.block(),block());
    }
    // assign with different layout (same dimRange)
    template <class Fy>
    void assign(const Derivatives<Fy,dimD,dimR,deriv,DerivativeLayoutNS::value> &y)
    {
      for (unsigned int rr=0; rr<dimR; ++rr)
        deriv_[rr].assign(y,rr);
    }
    // assign with scalar functions to component r
    template <class Fy,DerivativeLayoutNS::DerivativeLayout layouty>
    void assign(unsigned int r,const Derivatives<Fy,dimD,1,deriv,layouty> &y)
    {
      deriv_[r].assign(r,y);
    }

    Block &block()
    {
      return reinterpret_cast<Block&>(*this);
    }
    const Block &block() const
    {
      return reinterpret_cast<const Block&>(*this);
    }

    ScalarDeriv &operator[](int r) {
      return deriv_[r];
    }
    const ScalarDeriv &operator[](int r) const {
      return deriv_[r];
    }
  protected:
    Dune::FieldVector<ScalarDeriv,dimR> deriv_;
  };

  // ******************************************
  // AXPY *************************************
  // ******************************************
  template <class Vec1,class Vec2,unsigned int deriv>
  struct LFETensorAxpy
  {
    template <class Field>
    static void apply(unsigned int r,const Field &a,
                      const Vec1 &x, Vec2 &y)
    {
      y.axpy(a,x);
    }
  };
  template <class F1,int dimD,int dimR,
      unsigned int d,
      class Vec2,
      unsigned int deriv>
  struct LFETensorAxpy<Derivatives<F1,dimD,dimR,d,DerivativeLayoutNS::value>,Vec2,deriv>
  {
    typedef Derivatives<F1,dimD,dimR,d,DerivativeLayoutNS::value> Vec1;
    template <class Field>
    static void apply(unsigned int r,const Field &a,
                      const Vec1 &x, Vec2 &y)
    {
      const FieldVector<F1,Vec2::size> &xx = x.template block<deriv>();
      for (int i=0; i<y.size; ++i)
        y[i] += xx[i]*a;
    }
  };
  template <class F1,int dimD,int dimR,
      unsigned int d,
      class Vec2,
      unsigned int deriv>
  struct LFETensorAxpy<Derivatives<F1,dimD,dimR,d,DerivativeLayoutNS::derivative>,Vec2,deriv>
  {
    typedef Derivatives<F1,dimD,dimR,d,DerivativeLayoutNS::derivative> Vec1;
    template <class Field>
    static void apply(unsigned int r,const Field &a,
                      const Vec1 &x, Vec2 &y)
    {
      for (int rr=0; rr<dimR; ++rr)
        LFETensorAxpy<Derivatives<F1,dimD,1,d,DerivativeLayoutNS::value>,
            Vec2,deriv>::apply(rr,a,x[rr],y);
    }
  };
  template <class F1,int dimD,
      unsigned int d,
      class Vec2,
      unsigned int deriv>
  struct LFETensorAxpy<Derivatives<F1,dimD,1,d,DerivativeLayoutNS::derivative>,Vec2,deriv>
  {
    typedef Derivatives<F1,dimD,1,d,DerivativeLayoutNS::derivative> Vec1;
    template <class Field>
    static void apply(unsigned int r,const Field &a,
                      const Vec1 &x, Vec2 &y)
    {
      LFETensorAxpy<Derivatives<F1,dimD,1,d,DerivativeLayoutNS::value>,
          Vec2,deriv>::apply(r,a,x[0],y);
    }
  };
  template <class F1,int dimD,
      unsigned int d,
      class Vec2,
      unsigned int deriv>
  struct LFETensorAxpy<Derivatives<F1,dimD,1,d,DerivativeLayoutNS::value>,Vec2,deriv>
  {
    typedef Derivatives<F1,dimD,1,d,DerivativeLayoutNS::value> Vec1;
    template <class Field>
    static void apply(unsigned int r,const Field &a,
                      const Vec1 &x, Vec2 &y)
    {
      typedef LFETensor<F1,dimD,deriv> LFETensorType;
      const unsigned int rr = r*LFETensorType::size;
      const FieldVector<F1,LFETensorType::size> &xx = x.template block<deriv>();
      for (int i=0; i<FieldVector<F1,LFETensorType::size>::dimension; ++i)
        y[rr+i] += xx[i]*a;
    }
  };

  // ***********************************************
  // Assign ****************************************
  // ***********************************************
  template <class Vec1,class Vec2>
  struct DerivativeAssign
  {
    static void apply(unsigned int r,const Vec1 &vec1,Vec2 &vec2)
    {
      field_cast(vec1,vec2);
    }
  };
  template <int dimD,int dimR,unsigned int deriv, DerivativeLayoutNS::DerivativeLayout layout,
      class F1,class F2>
  struct DerivativeAssign<Derivatives<F1,dimD,dimR,deriv,layout>,
      Derivatives<F2,dimD,dimR,deriv,layout> >
  {
    typedef Derivatives<F1,dimD,dimR,deriv,layout> Vec1;
    typedef Derivatives<F2,dimD,dimR,deriv,layout> Vec2;
    static void apply(unsigned int r,const Vec1 &vec1,Vec2 &vec2)
    {
      field_cast(vec1.block(),vec2.block());
    }
  };
  template <int dimD,int dimR,unsigned int deriv,
      class F1, class F2>
  struct DerivativeAssign<Derivatives<F1,dimD,dimR,deriv,DerivativeLayoutNS::value>,
      Derivatives<F2,dimD,dimR,deriv,DerivativeLayoutNS::derivative> >
  {
    typedef Derivatives<F1,dimD,dimR,deriv,DerivativeLayoutNS::value> Vec1;
    typedef Derivatives<F2,dimD,dimR,deriv,DerivativeLayoutNS::derivative> Vec2;
    static void apply(unsigned int r,const Vec1 &vec1,Vec2 &vec2)
    {
      vec2.assign(vec1);
    }
  };
  template <int dimD,int dimR,unsigned int deriv,
      class F1, class F2>
  struct DerivativeAssign<Derivatives<F1,dimD,dimR,deriv,DerivativeLayoutNS::derivative>,
      Derivatives<F2,dimD,dimR,deriv,DerivativeLayoutNS::value> >
  {
    typedef Derivatives<F1,dimD,dimR,deriv,DerivativeLayoutNS::derivative> Vec1;
    typedef Derivatives<F2,dimD,dimR,deriv,DerivativeLayoutNS::value> Vec2;
    static void apply(unsigned int r,const Vec1 &vec1,Vec2 &vec2)
    {
      vec2.assign(vec1);
    }
  };
  template <int dimD,int dimR,unsigned int deriv,DerivativeLayoutNS::DerivativeLayout layout,
      class F1, class F2>
  struct DerivativeAssign<Derivatives<F1,dimD,1,deriv,layout>,
      Derivatives<F2,dimD,dimR,deriv,DerivativeLayoutNS::value> >
  {
    typedef Derivatives<F1,dimD,1,deriv,layout> Vec1;
    typedef Derivatives<F2,dimD,dimR,deriv,DerivativeLayoutNS::value> Vec2;
    static void apply(unsigned int r,const Vec1 &vec1,Vec2 &vec2)
    {
      vec2.assign(r,vec1);
    }
  };
  template <int dimD,int dimR,unsigned int deriv,DerivativeLayoutNS::DerivativeLayout layout,
      class F1, class F2>
  struct DerivativeAssign<Derivatives<F1,dimD,1,deriv,layout>,
      Derivatives<F2,dimD,dimR,deriv,DerivativeLayoutNS::derivative> >
  {
    typedef Derivatives<F1,dimD,1,deriv,layout> Vec1;
    typedef Derivatives<F2,dimD,dimR,deriv,DerivativeLayoutNS::derivative> Vec2;
    static void apply(unsigned int r,const Vec1 &vec1,Vec2 &vec2)
    {
      vec2.assign(r,vec1);
    }
  };
  template <int dimD,unsigned int deriv,
      class F1, class F2>
  struct DerivativeAssign<Derivatives<F1,dimD,1,deriv,DerivativeLayoutNS::value>,
      Derivatives<F2,dimD,1,deriv,DerivativeLayoutNS::value> >
  {
    typedef Derivatives<F1,dimD,1,deriv,DerivativeLayoutNS::value> Vec1;
    typedef Derivatives<F2,dimD,1,deriv,DerivativeLayoutNS::value> Vec2;
    static void apply(unsigned int r,const Vec1 &vec1,Vec2 &vec2)
    {
      field_cast(vec1.block(),vec2.block());
    }
  };
  template <int dimD,unsigned int deriv,
      class F1, class F2>
  struct DerivativeAssign<Derivatives<F1,dimD,1,deriv,DerivativeLayoutNS::derivative>,
      Derivatives<F2,dimD,1,deriv,DerivativeLayoutNS::derivative> >
  {
    typedef Derivatives<F1,dimD,1,deriv,DerivativeLayoutNS::derivative> Vec1;
    typedef Derivatives<F2,dimD,1,deriv,DerivativeLayoutNS::derivative> Vec2;
    static void apply(unsigned int r,const Vec1 &vec1,Vec2 &vec2)
    {
      field_cast(vec1.block(),vec2.block());
    }
  };
  template <int dimD,unsigned int deriv,
      class F1, class F2>
  struct DerivativeAssign<Derivatives<F1,dimD,1,deriv,DerivativeLayoutNS::derivative>,
      Derivatives<F2,dimD,1,deriv,DerivativeLayoutNS::value> >
  {
    typedef Derivatives<F1,dimD,1,deriv,DerivativeLayoutNS::derivative> Vec1;
    typedef Derivatives<F2,dimD,1,deriv,DerivativeLayoutNS::value> Vec2;
    static void apply(unsigned int r,const Vec1 &vec1,Vec2 &vec2)
    {
      field_cast(vec1.block(),vec2.block());
    }
  };
  template <int dimD,unsigned int deriv,
      class F1, class F2>
  struct DerivativeAssign<Derivatives<F1,dimD,1,deriv,DerivativeLayoutNS::value>,
      Derivatives<F2,dimD,1,deriv,DerivativeLayoutNS::derivative> >
  {
    typedef Derivatives<F1,dimD,1,deriv,DerivativeLayoutNS::value> Vec1;
    typedef Derivatives<F2,dimD,1,deriv,DerivativeLayoutNS::derivative> Vec2;
    static void apply(unsigned int r,const Vec1 &vec1,Vec2 &vec2)
    {
      field_cast(vec1.block(),vec2.block());
    }
  };
  template <int dimD,unsigned int deriv,DerivativeLayoutNS::DerivativeLayout layout,
      class F1, class F2>
  struct DerivativeAssign<Derivatives<F1,dimD,1,deriv,layout>,
      F2 >
  {
    typedef Derivatives<F1,dimD,1,deriv,layout> Vec1;
    typedef F2 Vec2;
    static void apply(unsigned int r,const Vec1 &vec1,Vec2 &vec2)
    {
      field_cast(vec1.block(),vec2);
    }
  };
  template <int dimD,int dimR,
      class F1,unsigned int deriv,
      class F2>
  struct DerivativeAssign<Derivatives<F1,dimD,dimR,deriv,DerivativeLayoutNS::value>,FieldVector<F2,dimR> >
  {
    typedef Derivatives<F1,dimD,dimR,deriv,DerivativeLayoutNS::value> Vec1;
    typedef FieldVector<F2,dimR> Vec2;
    static void apply(unsigned int r,const Vec1 &vec1,Vec2 &vec2)
    {
      field_cast(vec1.template block<0>(),vec2);
    }
  };
  template <int dimD,int dimR,
      class F1,unsigned int deriv,
      class F2>
  struct DerivativeAssign<Derivatives<F1,dimD,dimR,deriv,DerivativeLayoutNS::derivative>,FieldVector<F2,dimR> >
  {
    typedef Derivatives<F1,dimD,dimR,deriv,DerivativeLayoutNS::derivative> Vec1;
    typedef FieldVector<F2,dimR> Vec2;
    static void apply(unsigned int r,const Vec1 &vec1,Vec2 &vec2)
    {
      for (int rr=0; rr<dimR; ++rr)
        field_cast(vec1[rr].template tensor<0>()[0].block(),vec2[rr]);
    }
  };
  template <int dimD,
      class F1,unsigned int deriv,
      class F2,int dimR>
  struct DerivativeAssign<Derivatives<F1,dimD,1,deriv,DerivativeLayoutNS::value>,FieldVector<F2,dimR> >
  {
    typedef Derivatives<F1,dimD,1,deriv,DerivativeLayoutNS::value> Vec1;
    typedef FieldVector<F2,dimR> Vec2;
    static void apply(unsigned int r,const Vec1 &vec1,Vec2 &vec2)
    {
      field_cast(vec1.template tensor<0>()[0].block(),vec2[r]);
    }
  };
  template <int dimD,
      class F1,unsigned int deriv,
      class F2,int dimR>
  struct DerivativeAssign<Derivatives<F1,dimD,1,deriv,DerivativeLayoutNS::derivative>,FieldVector<F2,dimR> >
  {
    typedef Derivatives<F1,dimD,1,deriv,DerivativeLayoutNS::derivative> Vec1;
    typedef FieldVector<F2,dimR> Vec2;
    static void apply(unsigned int r,const Vec1 &vec1,Vec2 &vec2)
    {
      field_cast(vec1[0].template tensor<0>()[0].block(),vec2[r]);
    }
  };
  template <int dimD,
      class F1,unsigned int deriv,
      class F2>
  struct DerivativeAssign<Derivatives<F1,dimD,1,deriv,DerivativeLayoutNS::value>,FieldVector<F2,1> >
  {
    typedef Derivatives<F1,dimD,1,deriv,DerivativeLayoutNS::value> Vec1;
    typedef FieldVector<F2,1> Vec2;
    static void apply(unsigned int r,const Vec1 &vec1,Vec2 &vec2)
    {
      field_cast(vec1.template tensor<0>()[0].block(),vec2);
    }
  };
  template <int dimD,
      class F1,unsigned int deriv,
      class F2>
  struct DerivativeAssign<Derivatives<F1,dimD,1,deriv,DerivativeLayoutNS::derivative>,FieldVector<F2,1> >
  {
    typedef Derivatives<F1,dimD,1,deriv,DerivativeLayoutNS::derivative> Vec1;
    typedef FieldVector<F2,1> Vec2;
    static void apply(unsigned int /*r*/,const Vec1 &vec1,Vec2 &vec2)
    {
      field_cast(vec1[0].template tensor<0>()[0].block(),vec2);
    }
  };

  // ***********************************************
  // IO ********************************************
  // ***********************************************
  template <class F,int dimD,unsigned int deriv>
  std::ostream &operator<< ( std::ostream &out, const LFETensor< F,dimD,deriv > &tensor )
  {
    return out << tensor.block();
  }
#if 0
  template <class F,int dimD,unsigned int deriv>
  std::ostream &operator<< ( std::ostream &out, const ScalarDerivatives< F,dimD,deriv > &d )
  {
    out << static_cast<const ScalarDerivatives< F,dimD,deriv-1 > &>(d);
    out << " , " << d.tensor() << std::endl;
    return out;
  }
  template <class F,int dimD>
  std::ostream &operator<< ( std::ostream &out, const ScalarDerivatives< F,dimD,0 > &d )
  {
    out << d.tensor() << std::endl;
    return out;
  }
#endif
  template <class F,int dimD,int dimR,unsigned int deriv>
  std::ostream &operator<< ( std::ostream &out, const Derivatives< F,dimD,dimR,deriv,DerivativeLayoutNS::derivative > &d )
  {
    out << " ( ";
    out << d[0];
    for (int r=1; r<dimR; ++r)
    {
      out << " , " << d[r];
    }
    out << " ) " << std::endl;
    return out;
  }
  template <class F,int dimD,int dimR,unsigned int deriv>
  std::ostream &operator<< ( std::ostream &out, const Derivatives< F,dimD,dimR,deriv,DerivativeLayoutNS::value > &d )
  {
    out << static_cast<const Derivatives< F,dimD,dimR,deriv-1,DerivativeLayoutNS::value > &>(d);
    out << " ( ";
    out << d[0];
    for (int r=1; r<dimR; ++r)
    {
      out << " , " << d[r];
    }
    out << " ) " << std::endl;
    return out;
  }
  template <class F,int dimD,int dimR>
  std::ostream &operator<< ( std::ostream &out, const Derivatives< F,dimD,dimR,0,DerivativeLayoutNS::derivative > &d )
  {
    out << " ( ";
    out << d[0];
    for (int r=1; r<dimR; ++r)
    {
      out << " , " << d[r];
    }
    out << " ) " << std::endl;
    return out;
  }
  template <class F,int dimD,int dimR>
  std::ostream &operator<< ( std::ostream &out, const Derivatives< F,dimD,dimR,0,DerivativeLayoutNS::value > &d )
  {
    out << " ( ";
    out << d[0];
    for (int r=1; r<dimR; ++r)
    {
      out << " , " << d[r];
    }
    out << " ) " << std::endl;
    return out;
  }
  template <class F,int dimD,int dimR,unsigned int deriv,DerivativeLayoutNS::DerivativeLayout layout>
  std::ostream &operator<< ( std::ostream &out, const std::vector<Derivatives< F,dimD,dimR,deriv,layout > > &y )
  {
    out << "Number of basis functions: " << y.size() << std::endl;
    for (unsigned int i=0; i<y.size(); ++i)
    {
      out << "Base " << i << " : " << std::endl;
      out << y[i];
      out << std::endl;
    }
    return out;
  }
}
#endif // DUNE_TENSOR_HH
