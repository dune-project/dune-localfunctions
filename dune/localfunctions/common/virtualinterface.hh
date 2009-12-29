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

    //! @copydoc C0LocalBasisInterface::size
    unsigned int size () const
    {
      return impl_.size();
    }

    //! @copydoc C0LocalBasisInterface::order
    unsigned int order () const
    {
      return impl_.order();
    }

    //! @copydoc C0LocalBasisInterface::evaluateFunction
    inline void evaluateFunction (const typename Traits::DomainType& in,
                                  std::vector<typename Traits::RangeType>& out) const
    {
      impl_.evaluateFunction(in,out);
    }

  protected:
    const Imp impl_;

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

    typedef T Traits;
    using Base::size;
    using Base::order;
    using Base::evaluateFunction;

    //! @copydoc C0LocalBasisInterface::evaluateJacobian
    inline void evaluateJacobian(const typename Traits::DomainType& in,         // position
                                 std::vector<typename Traits::JacobianType>& out) const      // return value
    {
      impl_.evaluateJacobian(in,out);
    }

  protected:
    using Base::impl_;

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

    typedef typename Base::FunctionType FunctionType;

    typedef typename Base::CoefficientType CoefficientType;

    virtual void interpolate (const FunctionType& f, std::vector<CoefficientType>& out) const
    {
      impl_.interpolate(f,out);
    }

  protected:
    const Imp impl_;

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
  template<class LocalBasisTraits>
  class LocalFiniteElementVirtualInterface
  {
  public:
    //typedef T Traits;

    //! @copydoc LocalFiniteElementInterface::localBasis
    virtual const C0LocalBasisVirtualInterface<LocalBasisTraits>& localBasis () const = 0;

    //! @copydoc LocalFiniteElementInterface::localCoefficients
    virtual const LocalCoefficientsVirtualInterface& localCoefficients () const = 0;

    //! @copydoc LocalFiniteElementInterface::localInterpolation
    virtual const LocalInterpolationVirtualInterface<typename LocalBasisTraits::DomainType,
        typename LocalBasisTraits::RangeType>& localInterpolation () const = 0;

    //! @copydoc LocalFiniteElementInterface::type
    virtual const GeometryType type () const = 0;

  };


  /**
   * @brief class for wrapping a local finite element
   *        using the virtual interface
   *
   * @tparam Imp LocalFiniteElementInterface implementation
   **/
  template<class Imp>
  class LocalFiniteElementVirtualImp
    : public LocalFiniteElementVirtualInterface<typename Imp::Traits::LocalBasisType::Traits>
  {
    /** \brief The traits class of the local basis implementation.  Typedef'ed her for legibility */
    typedef typename Imp::Traits::LocalBasisType::Traits LocalBasisTraits;

  public:
    //typedef T Traits;

    //! @copydoc constructor taking a LocalFiniteElementInterface implementation
    LocalFiniteElementVirtualImp( const Imp &imp )
      : impl_(imp),
        localBasisImp_(impl_.localBasis()),
        localCoefficientsImp_(impl_.localCoefficients()),
        localInterpolationImp_(impl_.localInterpolation())
    {}

    //! Default constructor.  Assumes that the implementation class is default constructible as well.
    LocalFiniteElementVirtualImp()
      : impl_(),
        localBasisImp_(impl_.localBasis()),
        localCoefficientsImp_(impl_.localCoefficients()),
        localInterpolationImp_(impl_.localInterpolation())
    {}

    //! @copydoc LocalFiniteElementInterface::localBasis
    const C0LocalBasisVirtualInterface<LocalBasisTraits>& localBasis () const
    {
      return localBasisImp_;
    }

    //! @copydoc LocalFiniteElementInterface::localCoefficients
    const LocalCoefficientsVirtualInterface& localCoefficients () const
    {
      return localCoefficientsImp_;
    }

    //! @copydoc LocalFiniteElementInterface::localInterpolation
    const LocalInterpolationVirtualInterface<typename LocalBasisTraits::DomainType,
        typename LocalBasisTraits::RangeType>& localInterpolation () const
    {
      return localInterpolationImp_;
    }

    //! @copydoc LocalFiniteElementInterface::type
    const GeometryType type () const
    {
      return impl_.type();
    }

  protected:
    const Imp impl_;

    /** \todo This needs to automatically change to C0LocalBasisBla... to work with C0 shape functions */
    const C1LocalBasisVirtualImp<LocalBasisTraits, typename Imp::Traits::LocalBasisType> localBasisImp_;
    const LocalCoefficientsVirtualImp<typename Imp::Traits::LocalCoefficientsType> localCoefficientsImp_;
    const LocalInterpolationVirtualImp<typename LocalBasisTraits::DomainType,
        typename LocalBasisTraits::RangeType,
        typename Imp::Traits::LocalInterpolationType> localInterpolationImp_;


  };
}
#endif
