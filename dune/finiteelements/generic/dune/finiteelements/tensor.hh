// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
namespace Dune
{
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
      block = f;
      return *this;
    }
    This &operator=(const Block& b)
    {
      block = b;
      return *this;
    }
    Block block;
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
      block = f;
      return *this;
    }
    This &operator=(const Block& b)
    {
      block = b;
      return *this;
    }
    Block block;
  };
  // ****************************************
  // ****************************************
#if 0
  template <class F,int dimD,unsigned int deriv>
  struct ScalarDerivatives : public ScalarDerivatives<F,dimD,deriv-1>
  {
    typedef ScalarDerivatives<F,dimD,deriv-1> Base;

    template <unsigned int dorder>
    Tensor<F,dimD,dorder> &tensor()
    {
      return tensor(Int2Type<dorder>());
    }
  protected:
    template <unsigned int dorder>
    Tensor<F,dimD,dorder> &tensor(Int2Type<dorder> &dorderVar)
    {
      return Base::tensor(dorderVar);
    }
    Tensor<F,dimD,deriv> &tensor(Int2Type<deriv> &dorderVar)
    {
      return tensor_;
    }
    Tensor<F,dimD,deriv> tensor_;
  };
  template <class F,int dimD>
  struct ScalarDerivatives<F,dimD,0>
  {
    Tensor<F,dimD,0> &tensor()
    {
      return tensor_;
    }
    Tensor<F,dimD,0> tensor_;
  };
  // ***********************************************************
  template <class F,int dimD,int dimR,unsigned int deriv>
  struct Derivatives<F,dimD,deriv>
  {
    typedef ScalarDerivatives<F,dimD,deriv> ScalarDeriv;
    ScalarDeriv &operator[](int r) {
      return deriv_[r];
    }
  protected:
    Dune::FieldVector<ScalarDeriv,dimR> deriv_;
  };
#endif
  // *************************************************
  // *************************************************
  template <class F,int dimD,int dimR,unsigned int deriv>
  struct Derivative : public Derivative<F,dimD,dimR,deriv-1>
  {
    typedef Derivative<F,dimD,dimR,deriv> This;
    typedef Derivative<F,dimD,dimR,deriv-1> Base;
    typedef Tensor<F,dimD,deriv> ThisTensor;
    static const unsigned int size = Base::size+ThisTensor::size*dimR;

    This &operator=(const F& f)
    {
      Base::operator=(f);
      for (int r=0; r<dimR; ++r)
        tensor_[r] = f;
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
    template <int totalSize>
    This &operator=(const Dune::FieldVector<Dune::FieldVector<F,totalSize>,dimR> &t)
    {
      reinterpret_cast<Dune::FieldVector<Dune::FieldVector<F,totalSize>,dimR> &>(*this) = t;
      return *this;
    }
    template <int totalSize>
    This &operator=(const Dune::FieldVector<F,totalSize> &t)
    {
      reinterpret_cast<Dune::FieldVector<F,totalSize> &>(*this) = t;
      return *this;
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
    template <unsigned int dorder>
    Dune::FieldVector<Tensor<F,dimD,dorder>,dimR> &tensor(Int2Type<dorder> &dorderVar)
    {
      return Base::tensor(dorderVar);
    }
    Dune::FieldVector<Tensor<F,dimD,deriv>,dimR> &tensor(Int2Type<deriv> &dorderVar)
    {
      return tensor_;
    }
    Dune::FieldVector<ThisTensor,dimR> tensor_;
  };
  template <class F,int dimD,int dimR>
  struct Derivative<F,dimD,dimR,0>
  {
    typedef Derivative<F,dimD,dimR,0> This;
    typedef Tensor<F,dimD,0> ThisTensor;
    static const unsigned int size = ThisTensor::size*dimR;
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
    template <int totalSize>
    This &operator=(const Dune::FieldVector<Dune::FieldVector<F,totalSize>,dimR> &t)
    {
      reinterpret_cast<Dune::FieldVector<Dune::FieldVector<F,totalSize>,dimR> &>(*this) = t;
      return *this;
    }
    template <int totalSize>
    This &operator=(const Dune::FieldVector<F,totalSize> &t)
    {
      reinterpret_cast<Dune::FieldVector<F,totalSize> &>(*this) = t;
      return *this;
    }
    ThisTensor &operator[](int r) {
      return tensor_[r];
    }
    const ThisTensor &operator[](int r) const {
      return tensor_[r];
    }
    Dune::FieldVector<Tensor<F,dimD,0>,dimR> &tensor()
    {
      return tensor_;
    }
    Dune::FieldVector<ThisTensor,dimR> tensor_;
  };
  template <class F,int dimD,unsigned int deriv>
  std::ostream &operator<< ( std::ostream &out, const Tensor< F,dimD,deriv > &tensor )
  {
    return out << tensor.block;
  }
  template <class F,int dimD,int dimR,unsigned int deriv>
  std::ostream &operator<< ( std::ostream &out, const Derivative< F,dimD,dimR,deriv > &d )
  {
    out << Derivative< F,dimD,dimR,deriv-1 >(d);
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
  std::ostream &operator<< ( std::ostream &out, const Derivative< F,dimD,dimR,0 > &d )
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
  std::ostream &operator<< ( std::ostream &out, const std::vector<Derivative< F,dimD,dimR,deriv > > &y )
  {
    out << "Number of basis functions: " << y.size() << std::endl;
    for (int i=0; i<y.size(); ++i)
    {
      out << "Base " << i << " : " << std::endl;
      out << y[i];
      out << std::endl;
    }
    return out;
  }
}
