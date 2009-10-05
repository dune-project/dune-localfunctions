// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_TENSOR_HH
#define DUNE_TENSOR_HH

#include <dune/common/fvector.hh>
#include <dune/common/misc.hh>

namespace Dune
{
  // Structure for scalar tensor of order deriv
  template <class F,int dimD,unsigned int deriv>
  struct Tensor
  {
    typedef Tensor<F,dimD,deriv> This;
    typedef Tensor<F,dimD-1,deriv> BaseDim;
    typedef Tensor<F,dimD,deriv-1> BaseDeriv;
    static const unsigned int size = BaseDim::size+BaseDeriv::size;
    typedef Dune::FieldVector<F,size> Block;
    This &operator=(const F& f)
    {
      block() = f;
      return *this;
    }
    This &operator=(const Block& b)
    {
      block() = b;
      return *this;
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
    void assign(const Tensor<Fy,dimD,deriv> &y)
    {
      field_cast(y.block(),block());
    }
    Block block_;
  };
  // ******************************************
  template <class F,unsigned int deriv>
  struct Tensor<F,0,deriv>
  {
    static const int size = 0;
  };
  template <class F>
  struct Tensor<F,0,0>
  {
    static const int size = 1;
  };
  template <class F,int dimD>
  struct Tensor<F,dimD,0>
  {
    typedef Tensor<F,dimD,0> This;
    static const int size = 1;
    typedef Dune::FieldVector<F,size> Block;
    This &operator=(const F& f)
    {
      block() = f;
      return *this;
    }
    This &operator=(const Block& b)
    {
      block() = b;
      return *this;
    }
    void axpy(const F& a, const This &y)
    {
      block().axpy(a,y.block());
    }
    template <class Fy>
    void assign(const Tensor<Fy,dimD,0> &y)
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
  enum DerivativeLayout {value,derivative};
  template <class F,int dimD,int dimR,unsigned int deriv,
      DerivativeLayout layout>
  struct Derivatives;

  // Implemnetation for valued based layout
  template <class F,int dimD,int dimR,unsigned int deriv>
  struct Derivatives<F,dimD,dimR,deriv,value>
    : public Derivatives<F,dimD,dimR,deriv-1,value>
  {
    typedef Derivatives<F,dimD,dimR,deriv,value> This;
    typedef Derivatives<F,dimD,dimR,deriv-1,value> Base;
    typedef Tensor<F,dimD,deriv> ThisTensor;
    static const DerivativeLayout layout = value;
    static const unsigned int dimDomain = dimD;
    static const unsigned int dimRange = dimR;
    static const unsigned int size = Base::size+ThisTensor::size*dimR;
    typedef Dune::FieldVector<F,size> Block;

    This &operator=(const F& f)
    {
      block() = f;
      return *this;
    }
    This &operator=(const Dune::FieldVector<ThisTensor,dimR> &t)
    {
      tensor_ = t;
      return *this;
    }
    template <unsigned int dorder>
    This &operator=(const Dune::FieldVector<Tensor<F,dimD,dorder>,dimR> &t)
    {
      tensor<dorder>() = t;
      return *this;
    }
    This &operator=(const Block &t)
    {
      block() = t;
      return *this;
    }

    void axpy(const F &a, const This &y)
    {
      block().axpy(a,y.block());
    }
    // assign with same layout (only diffrent Field)
    template <class Fy>
    void assign(unsigned int r,const Derivatives<Fy,dimD,dimR,deriv,value> &y)
    {
      field_cast(y.block(),block());
    }
    // assign with diffrent layout (same dimRange)
    template <class Fy>
    void assign(unsigned int r,const Derivatives<Fy,dimD,dimR,deriv,derivative> &y)
    {
      Base::assign(r,y);
      for (int rr=0; rr<dimR; ++rr)
        tensor_[rr] = y[rr].template tensor<deriv>()[0];
    }
    // assign with scalar functions to component r
    template <class Fy,DerivativeLayout layouty>
    void assign(unsigned int r,const Derivatives<Fy,dimD,1,deriv,layouty> &y)
    {
      Base::assign(r,static_cast<const Derivatives<Fy,dimD,1,deriv-1,layouty>&>(y));
      tensor_[r].assign(y[0]);
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
    const Dune::FieldVector<Tensor<F,dimD,dorder>,dimR> &tensor() const
    {
      const Int2Type<dorder> a = Int2Type<dorder>();
      return tensor(a);
    }
    template <unsigned int dorder>
    Dune::FieldVector<Tensor<F,dimD,dorder>,dimR> &tensor()
    {
      return tensor(Int2Type<dorder>());
    }
    ThisTensor &operator[](int r) {
      return tensor_[r];
    }
    const ThisTensor &operator[](int r) const {
      return tensor_[r];
    }
  protected:
    template <class Fy,unsigned int dy>
    void assign(unsigned int r,const Derivatives<Fy,dimD,dimR,dy,derivative> &y)
    {
      Base::assign(r,y);
      for (int rr=0; rr<dimR; ++rr)
        tensor_[rr] = y[rr].template tensor<deriv>()[0];
    }
    template <class Fy,unsigned int dy,DerivativeLayout layouty>
    void assign(unsigned int r,const Derivatives<Fy,dimD,1,dy,layouty> &y)
    {
      Base::assign(r,y);
      tensor_[r].assign(y[r]);
    }
    template <int dorder>
    const Dune::FieldVector<Tensor<F,dimD,dorder>,dimR> &tensor(const Int2Type<dorder> &dorderVar) const
    {
      return Base::tensor(dorderVar);
    }
    const Dune::FieldVector<Tensor<F,dimD,deriv>,dimR> &tensor(const Int2Type<deriv> &dorderVar) const
    {
      return tensor_;
    }
    template <int dorder>
    Dune::FieldVector<Tensor<F,dimD,dorder>,dimR> &tensor(const Int2Type<dorder> &dorderVar)
    {
      return Base::tensor(dorderVar);
    }
    Dune::FieldVector<Tensor<F,dimD,deriv>,dimR> &tensor(const Int2Type<deriv> &dorderVar)
    {
      return tensor_;
    }
    Dune::FieldVector<ThisTensor,dimR> tensor_;
  };
  template <class F,int dimD,int dimR>
  struct Derivatives<F,dimD,dimR,0,value>
  {
    typedef Derivatives<F,dimD,dimR,0,value> This;
    typedef Tensor<F,dimD,0> ThisTensor;
    static const DerivativeLayout layout = value;
    static const unsigned int dimDomain = dimD;
    static const unsigned int dimRange = dimR;
    static const unsigned int size = ThisTensor::size*dimR;
    typedef Dune::FieldVector<F,size> Block;

    This &operator=(const F& f)
    {
      for (int r=0; r<dimR; ++r)
        tensor_[r] = f;
      return *this;
    }
    This &operator=(const Dune::FieldVector<ThisTensor,dimR> &t)
    {
      tensor_ = t;
      return *this;
    }

    This &operator=(const Block &t)
    {
      block() = t;
      return *this;
    }
    void axpy(const F &a, const This &y)
    {
      block().axpy(a,y.block());
    }
    template <class Fy>
    void assign(unsigned int r,const Derivatives<Fy,dimD,dimR,0,value> &y)
    {
      field_cast(y.block(),block());
    }
    template <class Fy>
    void assign(unsigned int r,const Derivatives<Fy,dimD,dimR,0,derivative> &y)
    {
      for (int rr=0; rr<dimR; ++rr)
        tensor_[rr] = y[rr].template tensor<0>()[0];
    }
    template <class Fy,DerivativeLayout layouty>
    void assign(unsigned int r,const Derivatives<Fy,dimD,1,0,layouty> &y)
    {
      tensor_[r].assign(y[r]);
    }

    Block &block()
    {
      return reinterpret_cast<Block&>(*this);
    }
    const Block &block() const
    {
      return reinterpret_cast<const Block&>(*this);
    }

    ThisTensor &operator[](int r) {
      return tensor_[r];
    }
    const ThisTensor &operator[](int r) const {
      return tensor_[r];
    }
    template <int dorder>
    const Dune::FieldVector<Tensor<F,dimD,0>,dimR> &tensor() const
    {
      return tensor_;
    }
    Dune::FieldVector<Tensor<F,dimD,0>,dimR> &tensor()
    {
      return tensor_;
    }
  protected:
    const Dune::FieldVector<Tensor<F,dimD,0>,dimR> &tensor(const Int2Type<0> &dorderVar) const
    {
      return tensor_;
    }
    Dune::FieldVector<Tensor<F,dimD,0>,dimR> &tensor(const Int2Type<0> &dorderVar)
    {
      return tensor_;
    }
    template <class Fy,unsigned int dy>
    void assign(unsigned int r,const Derivatives<Fy,dimD,dimR,dy,derivative> &y)
    {
      for (int rr=0; rr<dimR; ++rr)
        tensor_[rr] = y[rr].template tensor<0>()[0];
    }
    template <class Fy,unsigned int dy,DerivativeLayout layouty>
    void assign(unsigned int r,const Derivatives<Fy,dimD,1,dy,layouty> &y)
    {
      tensor_[r].assign(y[r]);
    }
    Dune::FieldVector<ThisTensor,dimR> tensor_;
  };
  // Implemnetation for derivative based layout
  template <class F,int dimD,int dimR,unsigned int deriv>
  struct Derivatives<F,dimD,dimR,deriv,derivative>
  {
    typedef Derivatives<F,dimD,dimR,deriv,derivative> This;
    typedef Derivatives<F,dimD,1,deriv,value> ScalarDeriv;
    static const DerivativeLayout layout = value;
    static const unsigned int dimDomain = dimD;
    static const unsigned int dimRange = dimR;
    static const unsigned int size = ScalarDeriv::size*dimR;
    typedef Dune::FieldVector<F,size> Block;

    This &operator=(const F& f)
    {
      block() = f;
      return *this;
    }
    This &operator=(const Block &t)
    {
      block() = t;
      return *this;
    }

    void axpy(const F &a, const This &y)
    {
      block().axpy(a,y.block());
    }
    // assign with same layout (only diffrent Field)
    template <class Fy>
    void assign(unsigned int r,const Derivatives<Fy,dimD,dimR,deriv,derivative> &y)
    {
      field_cast(y.block(),block());
    }
    // assign with diffrent layout (same dimRange)
    template <class Fy>
    void assign(unsigned int r,const Derivatives<Fy,dimD,dimR,deriv,value> &y)
    {
      abort();
    }
    // assign with scalar functions to component r
    template <class Fy,DerivativeLayout layouty>
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

  // ***********************************************
  // IO *********************************************
  // ***********************************************
  template <class F,int dimD,unsigned int deriv>
  std::ostream &operator<< ( std::ostream &out, const Tensor< F,dimD,deriv > &tensor )
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
  std::ostream &operator<< ( std::ostream &out, const Derivatives< F,dimD,dimR,deriv,derivative > &d )
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
  std::ostream &operator<< ( std::ostream &out, const Derivatives< F,dimD,dimR,deriv,value > &d )
  {
    out << static_cast<const Derivatives< F,dimD,dimR,deriv-1,value > &>(d);
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
  std::ostream &operator<< ( std::ostream &out, const Derivatives< F,dimD,dimR,0,derivative > &d )
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
  std::ostream &operator<< ( std::ostream &out, const Derivatives< F,dimD,dimR,0,value > &d )
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
  template <class F,int dimD,int dimR,unsigned int deriv,DerivativeLayout layout>
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
