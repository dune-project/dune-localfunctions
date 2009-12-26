// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_VIRTUALINTERFACE_HH
#define DUNE_VIRTUALINTERFACE_HH

#include <dune/common/function.hh>

#include <dune/localfunctions/common/localbasis.hh>
#include <dune/localfunctions/common/localinterpolation.hh>
#include <dune/localfunctions/common/localcoefficients.hh>
#include <dune/localfunctions/common/localfiniteelement.hh>

namespace Dune
{
  // -----------------------------------------------------------------
  // Basis
  // -----------------------------------------------------------------
  /**
   * @brief virtual base class for a C0 local basis
   *
   * This class defines the same interface as the C0LocalBasisInterface
   * class but using pure virtual methods
   **/
  template<class T>
  class C0LocalBasisVirtualInterface
  {
  public:
    typedef T Traits;

    //! @copydoc C0LocalBasisInterface::size
    virtual unsigned int size () const = 0;

    //! @copydoc C0LocalBasisInterface::order
    virtual unsigned int order () const = 0;

    //! @copydoc C0LocalBasisInterface::evaluateFunction
    virtual void evaluateFunction (const typename Traits::DomainType& in,
                                   std::vector<typename Traits::RangeType>& out) const = 0;
  };


  /**
   * @brief class for wrapping a C0 basis using the virtual interface
   *
   * @tparam T C0LocalBasisTraits class
   * @tparam Imp C0LocalBasisInterface implementation
   **/
  template<class T , class Imp>
  class C0LocalBasisVirtualImp
    : public virtual C0LocalBasisVirtualInterface<T>
  {
  public:
    typedef T Traits;

    //! constructor taking a Dune::C0LocalBasisInterface implementation
    C0LocalBasisVirtualImp( const Imp &imp )
      : impl_(imp) {}

    //! constructor taking a Dune::LocalFiniteElementInterface implementation
    //! the local basis is taken from the finite element
    template <class FET,class FEImpl>
    C0LocalBasisVirtualImp( const LocalFiniteElementInterface<FET, FEImpl> &fe )
      : impl_(fe.localBasis()) {}

    //! @copydoc C0LocalBasisInterface::size
    unsigned int size () const
    {
      return impl().size();
    }

    //! @copydoc C0LocalBasisInterface::order
    unsigned int order () const
    {
      return impl().order();
    }

    //! @copydoc C0LocalBasisInterface::evaluateFunction
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

  /**
   * @brief virtual base class for a C1 local basis
   *
   * This class defines the same interface as the C1LocalBasisInterface
   * class but using pure virtual methods
   **/
  template<class T>
  class C1LocalBasisVirtualInterface
    : public virtual C0LocalBasisVirtualInterface<T>
  {
  public:

    typedef T Traits;
    typedef C0LocalBasisVirtualInterface<T> Base;
    using Base::size;
    using Base::order;
    using Base::evaluateFunction;

    //! @copydoc C0LocalBasisInterface::evaluateJacobian
    virtual void evaluateJacobian(const typename Traits::DomainType& in,         // position
                                  std::vector<typename Traits::JacobianType>& out) const = 0;
  };
  /**
   * @brief class for wrapping a C1 basis using the virtual interface
   *
   * @tparam T C1LocalBasisTraits class
   * @tparam Imp C1LocalBasisInterface implementation
   **/
  template<class T , class Imp>
  class C1LocalBasisVirtualImp
    : public virtual C1LocalBasisVirtualInterface<T>
      , public virtual C0LocalBasisVirtualImp<T,Imp>
  {
    typedef C0LocalBasisVirtualImp<T,Imp> Base;

  public:
    //! constructor taking a Dune::C1LocalBasisInterface implementation
    C1LocalBasisVirtualImp( const Imp &imp )
      : Base(imp)
    {}

    //! constructor taking a Dune::LocalFiniteElementInterface implementation
    //! the local basis is taken from the finite element
    template <class FET,class FEImpl>
    C1LocalBasisVirtualImp( const LocalFiniteElementInterface<FET, FEImpl> &fe )
      : Base(fe.localBasis())
    {}

    typedef T Traits;
    using Base::size;
    using Base::order;
    using Base::evaluateFunction;

    //! @copydoc C0LocalBasisInterface::evaluateJacobian
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
  /**
   * @brief virtual base class for a local interpolation
   *
   * This class defines the same interface as the LocalInterpolationInterface
   * class but using pure virtual methods
   **/
  template<class DomainType, class RangeType>
  class LocalInterpolationVirtualInterface
  {
  public:

    //! type of virtual function to interpolate
    typedef Dune::VirtualFunction<DomainType, RangeType> FunctionType;

    //! type of the coefficient vector in the interpolate method
    typedef typename RangeType::field_type CoefficientType;

    //! @copydoc LocalInterpolationInterface::interpolate
    //! this is the pure virtual method taking a VirtualFunction
    virtual void interpolate (const FunctionType& f, std::vector<CoefficientType>& out) const = 0;

    //! @copydoc LocalInterpolationInterface::interpolate
    //! this uses the pure virtual method by wrapping the template argument into a VirtualFunction
    template<typename F>
    void interpolate (const F& f, std::vector<CoefficientType>& out) const
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


  /**
   * @brief class for wrapping a local interpolation
   *        using the virtual interface
   *
   * @tparam DomainType domain type of the Dune::VirtualFunction to interpolate
   * @tparam RangeType range type of the Dune::VirtualFunction to interpolate
   * @tparam Imp LocalInterpolationInterface implementation
   **/
  template<class DomainType, class RangeType, class Imp>
  class LocalInterpolationVirtualImp
    : public LocalInterpolationVirtualInterface< DomainType, RangeType >
  {
    typedef LocalInterpolationVirtualInterface< DomainType, RangeType > Base;

  public:

    //! constructor taking a Dune::LocalInterpolationInterface implementation
    LocalInterpolationVirtualImp( const Imp &imp)
      : impl_(imp) {}

    //! constructor taking a Dune::LocalFiniteElementInterface implementation
    //! the local interpolation is taken from the finite element
    template <class FET,class FEImpl>
    LocalInterpolationVirtualImp( const LocalFiniteElementInterface<FET, FEImpl> &fe )
      : impl_(fe.localInterpolation()) {}

    typedef typename Base::FunctionType FunctionType;

    typedef typename Base::CoefficientType CoefficientType;

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
  /**
   * @brief virtual base class for local coefficients
   *
   * This class defines the same interface as the LocalCoefficientsInterface
   * class but using pure virtual methods
   **/
  class LocalCoefficientsVirtualInterface
  {
  public:

    //! @copydoc LocalCoefficientsInterface::size
    virtual std::size_t size () const = 0;

    //! @copydoc LocalCoefficientsInterface::localKey
    const virtual LocalKey& localKey (std::size_t i) const = 0;

  };


  /**
   * @brief class for wrapping local coefficients
   *        using the virtual interface
   *
   * @tparam Imp LocalCoefficientsInterface implementation
   **/
  template<class Imp>
  class LocalCoefficientsVirtualImp
    : public LocalCoefficientsVirtualInterface
  {
  public:
    //! constructor taking a Dune::LocalCoefficientsInterface implementation
    LocalCoefficientsVirtualImp( const Imp &imp )
      : impl_(imp) {}

    //! constructor taking a Dune::LocalFiniteElementInterface implementation
    //! the local coefficients is taken from the finite element
    template <class FET,class FEImpl>
    LocalCoefficientsVirtualImp( const LocalFiniteElementInterface<FET, FEImpl> &fe )
      : impl_(fe.localCoefficients()) {}

    //! @copydoc LocalCoefficientsInterface::size
    std::size_t size () const
    {
      return impl_.size();
    }

    //! @copydoc LocalCoefficientsInterface::localKey
    const LocalKey& localKey (std::size_t i) const
    {
      return impl_.localKey(i);
    }

  protected:
    const Imp impl_;

  };


  // -----------------------------------------------------------------
  // Finite Element
  // -----------------------------------------------------------------
  /**
   * @brief virtual base class for local finite elements
   *
   * This class defines the same interface as the LocalFiniteElementInterface
   * class but using pure virtual methods
   **/
  template<class T>
  class LocalFiniteElementVirtualInterface
  {
  public:
    typedef T Traits;

    //! @copydoc LocalFiniteElementInterface::localBasis
    virtual const typename T::LocalBasisType& localBasis () const = 0;

    //! @copydoc LocalFiniteElementInterface::localCoefficients
    virtual const LocalCoefficientsVirtualInterface& localCoefficients () const = 0;

    //! @copydoc LocalFiniteElementInterface::localInterpolation
    virtual const typename T::LocalInterpolationType& localInterpolation () const = 0;

    //! @copydoc LocalFiniteElementInterface::type
    virtual const GeometryType type () const = 0;

  };


  /**
   * @brief class for wrapping a local finite element
   *        using the virtual interface
   *
   * @tparam T traits class for the local finite element
   * @tparam Imp LocalFiniteElementInterface implementation
   **/
  template<class T, class Imp>
  class LocalFiniteElementVirtualImp
    : public LocalFiniteElementVirtualInterface<T>
  {
  public:
    typedef T Traits;

    //! @copydoc constructor taking a LocalFiniteElementInterface implementation
    LocalFiniteElementVirtualImp( const Imp &imp )
      : impl_(imp) {}

    //! Default constructor.  Assumes that the implementation class is default constructible as well.
    LocalFiniteElementVirtualImp()
    {}

    //! @copydoc LocalFiniteElementInterface::localBasis
    const typename T::LocalBasisType& localBasis () const
    {
      return impl_.localBasis();
    }

    //! @copydoc LocalFiniteElementInterface::localCoefficients
    const typename T::LocalCoefficientsType& localCoefficients () const
    {
      return impl_.localCoefficients();
    }

    //! @copydoc LocalFiniteElementInterface::localInterpolation
    const typename T::LocalInterpolationType& localInterpolation () const
    {
      return impl_.localInterpolation();
    }

    //! @copydoc LocalFiniteElementInterface::type
    const GeometryType type () const
    {
      return impl_.type();
    }

  protected:
    const Imp impl_;

  };
}
#endif
