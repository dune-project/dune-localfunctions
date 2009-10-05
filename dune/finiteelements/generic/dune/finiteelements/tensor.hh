// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
namespace Dune
{
  template <class F,int dimD,unsigned int deriv>
  struct Tensor
  {
    typedef Tensor<F,dimD-1,deriv> BaseDim;
    typedef Tensor<F,dimD,deriv-1> BaseDeriv;
    static const int size = BaseDim::size+BaseDeriv;
    typedef Dune::FieldVector<F,size> Block;
    Block block;
  };
  // ******************************************
  template <class F,unsigned int deriv>
  struct Tensor<F,0,deriv>
  {
    static const int size = 0;
  };
  template <class F,unsigned int deriv>
  struct Tensor<F,0,0>
  {
    static const int size = 1;
  };
  template <class F,int dimD>
  struct Tensor<F,dim,0>
  {
    static const int size = 1;
    typedef Dune::FieldVector<F,size> Block;
    Block block;
  };
  // ****************************************
  // ****************************************
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
  // *************************************************
  // *************************************************
  template <class F,int dimD,int dimR,unsigned int deriv>
  struct Derivative : public Derivative<F,dimD,dimR,deriv-1>
  {
    typedef Derivative<F,dimD,dimR,deriv-1> Base;

    template <unsigned int dorder>
    Dune::FieldVector<Tensor<F,dimD,dorder>,dimR> &tensor()
    {
      return tensor(Int2Type<dorder>());
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
    Dune::FieldVector<Tensor<F,dimD,deriv>,dimR> tensor_;
  };
  template <class F,int dimD,int dimR>
  struct Derivative<F,dimD,dimR,0>
  {
    Dune::FieldVector<Tensor<F,dimD,0>,dimR> &tensor()
    {
      return tensor_;
    }
    Dune::FieldVector<Tensor<F,dimD,0>,dimR> tensor_;
  };
}
