// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_VIRTUALINTERFACE_HH
#define DUNE_VIRTUALINTERFACE_HH

#include <dune/common/function.hh>

#include <dune/finiteelements/common/localbasis.hh>
#include <dune/finiteelements/common/localinterpolation.hh>
#include <dune/finiteelements/common/localcoefficients.hh>
#include <dune/finiteelements/common/localfiniteelement.hh>

namespace Dune
{
  // -----------------------------------------------------------------
  // Basis
  // -----------------------------------------------------------------
  // C0
  template<class T>
  class C0LocalBasisVirtualInterface
    : public C0LocalBasisInterface< T, C0LocalBasisVirtualInterface<T> >
  {
  public:
    typedef T Traits;
    virtual unsigned int size () const = 0;
    virtual unsigned int order () const = 0;
    virtual void evaluateFunction (const typename Traits::DomainType& in,
                                   std::vector<typename Traits::RangeType>& out) const = 0;
  };
  template<class T , class Imp>
  class C0LocalBasisVirtualImp
    : public virtual C0LocalBasisVirtualInterface<T>
  {
  public:
    typedef T Traits;

    C0LocalBasisVirtualImp( const Imp &imp )
      : impl_(imp) {}
    template <class FET,class FEImpl>
    C0LocalBasisVirtualImp( const LocalFiniteElementInterface<FET, FEImpl> &fe )
      : impl_(fe.localBasis()) {}
    unsigned int size () const
    {
      return impl().size();
    }
    unsigned int order () const
    {
      return impl().order();
    }
    inline void evaluateFunction (const typename Traits::DomainType& in,
                                  std::vector<typename Traits::RangeType>& out) const
    {
      impl().evaluateFunction(in,out);
    }
  protected:
    const Imp impl_;
  private:
    typedef C0LocalBasisInterface<T,Imp> Interface;
    const Interface& impl () const {return impl_;}
  };

  // C1
  template<class T>
  class C1LocalBasisVirtualInterface
    : public C1LocalBasisInterface< T, C1LocalBasisVirtualInterface<T> >
      , public virtual C0LocalBasisVirtualInterface<T>
  {
  public:
    typedef T Traits;
    typedef C0LocalBasisVirtualInterface<T> Base;
    using Base::size;
    using Base::order;
    using Base::evaluateFunction;
    virtual void evaluateJacobian(const typename Traits::DomainType& in,         // position
                                  std::vector<typename Traits::JacobianType>& out) const = 0;
  };
  template<class T , class Imp>
  class C1LocalBasisVirtualImp
    : public virtual C1LocalBasisVirtualInterface<T>
      , public virtual C0LocalBasisVirtualImp<T,Imp>
  {
    typedef C0LocalBasisVirtualImp<T,Imp> Base;
  public:
    C1LocalBasisVirtualImp( const Imp &imp )
      : Base(imp)
    {}
    template <class FET,class FEImpl>
    C1LocalBasisVirtualImp( const LocalFiniteElementInterface<FET, FEImpl> &fe )
      : Base(fe.localBasis())
    {}
    typedef T Traits;
    using Base::size;
    using Base::order;
    using Base::evaluateFunction;
    inline void evaluateJacobian(const typename Traits::DomainType& in,         // position
                                 std::vector<typename Traits::JacobianType>& out) const      // return value
    {
      impl().evaluateJacobian(in,out);
    }
  protected:
    using Base::impl_;
  private:
    typedef C1LocalBasisInterface<T,Imp> Interface;
    const Interface& impl () const {return impl_;}
  };

  // -----------------------------------------------------------------
  // Interpolation
  // -----------------------------------------------------------------
  template<class DomainType, class RangeType>
  class LocalInterpolationVirtualInterface
    : public LocalInterpolationInterface< LocalInterpolationVirtualInterface<DomainType,RangeType> >
  {
  public:
    typedef Dune::VirtualFunctionInterface<DomainType, RangeType> FunctionType;
    typedef typename RangeType::field_type CoefficientType;

    virtual void interpolate (const FunctionType& f, std::vector<CoefficientType>& out) const = 0;
    template<typename F, typename C> // is this C a good idea - in the basis function this type comes from the traits class
    void interpolate (const F& f, std::vector<C>& out) const
    {
      interpolate(VirtualFunctionWrapper<F>(f),out);
    }

  private:
    template <typename F>
    struct VirtualFunctionWrapper
      : public FunctionType
    {
    public:
      VirtualFunctionWrapper(const F &f) : f_(f) {}
      virtual void evaluate(const DomainType& x, RangeType& y) const
      {
        f_.evaluate(x,y);
      }
      const F &f_;
    };
  };
  template<class DomainType, class RangeType, class Imp>
  class LocalInterpolationVirtualImp
    : public LocalInterpolationVirtualInterface< DomainType, RangeType >
  {
  public:
    LocalInterpolationVirtualImp( const Imp &imp)
      : impl_(imp) {}
    template <class FET,class FEImpl>
    LocalInterpolationVirtualImp( const LocalFiniteElementInterface<FET, FEImpl> &fe )
      : impl_(fe.localInterpolation()) {}
    typedef VirtualFunctionInterface<DomainType, RangeType> FunctionType;
    typedef typename RangeType::field_type CoefficientType;

    virtual void interpolate (const FunctionType& f, std::vector<CoefficientType>& out) const
    {
      impl().interpolate(f,out);
    }
  protected:
    const Imp impl_;
  private:
    typedef LocalInterpolationInterface<Imp> Interface;
    const Interface& impl () const {return impl_;}
  };

  // -----------------------------------------------------------------
  // Coefficients
  // -----------------------------------------------------------------
  class LocalCoefficientsVirtualInterface
    : public LocalCoefficientsInterface< LocalCoefficientsVirtualInterface >
  {
  public:
    virtual std::size_t size () const = 0;
    const virtual LocalKey& localKey (std::size_t i) const = 0;
  };
  template<class Imp>
  class LocalCoefficientsVirtualImp
    : public LocalCoefficientsVirtualInterface
  {
  public:
    LocalCoefficientsVirtualImp( const Imp &imp )
      : impl_(imp) {}
    template <class FET,class FEImpl>
    LocalCoefficientsVirtualImp( const LocalFiniteElementInterface<FET, FEImpl> &fe )
      : impl_(fe.localCoefficients()) {}
    std::size_t size () const
    {
      return impl().size();
    }
    const LocalKey& localKey (std::size_t i) const
    {
      return impl().localKey(i);
    }
  protected:
    const Imp impl_;
  private:
    typedef LocalCoefficientsInterface<Imp> Interface;
    const Interface& impl () const {return impl_;}
  };

  // -----------------------------------------------------------------
  // Finite Element
  // -----------------------------------------------------------------
  template<class T>
  class LocalFiniteElementVirtualInterface
    : public LocalFiniteElementInterface< T, LocalFiniteElementVirtualInterface<T> >
  {
  public:
    typedef T Traits;
    LocalFiniteElementVirtualInterface( const GeometryType &gt ) : gt_(gt) {}
    virtual const typename T::LocalBasisType& localBasis () const = 0;
    virtual const typename T::LocalCoefficientsType& localCoefficients () const = 0;
    virtual const typename T::LocalInterpolationType& localInterpolation () const = 0;
    const GeometryType &type () const
    {
      return gt_;
    }
  private:
    const GeometryType gt_;
  };
  template<class T, class Imp>
  class LocalFiniteElementVirtualImp
    : public LocalFiniteElementVirtualInterface<T>
  {
  public:
    typedef T Traits;
    LocalFiniteElementVirtualImp( const Imp &imp )
      : impl_(imp) {}
    const typename T::LocalBasisType& localBasis () const
    {
      return impl().localBasis();
    }
    const typename T::LocalCoefficientsType& localCoefficients () const
    {
      return impl().localCoefficients();
    }
    const typename T::LocalInterpolationType& localInterpolation () const
    {
      return impl().localInterpolation();
    }
  protected:
    const Imp impl_;
  private:
    typedef LocalFiniteElementInterface<T,Imp> Interface;
    const Interface& impl () const {return impl_;}
  };
}
#endif
