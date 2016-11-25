// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_VIRTUALWRAPPERS_HH
#define DUNE_VIRTUALWRAPPERS_HH

#include <array>

#include <dune/common/deprecated.hh>
#include <dune/common/function.hh>

#include <dune/localfunctions/common/localbasis.hh>
#include <dune/localfunctions/common/localkey.hh>
#include <dune/localfunctions/common/virtualinterface.hh>

namespace Dune
{

  // forward declaration needed by friend declarations
  template<class Imp>
  class LocalFiniteElementVirtualImp;

  // default clone method is the copy constructor
  template<class Imp, bool IsInterface>
  struct LocalFiniteElementCloneFactoryHelper
  {
    static Imp* clone(const Imp& imp)
    {
      return new Imp(imp);
    }
  };

  // if FE derives from virtual interface the clone method is used
  template<class Imp>
  struct LocalFiniteElementCloneFactoryHelper<Imp, true>
  {
    static Imp* clone(const Imp& imp)
    {
      return imp.clone();
    }
  };

  // factory template to clone and create an objects
  template<class Imp>
  struct LocalFiniteElementCloneFactory
  {
    typedef LocalFiniteElementVirtualInterface<typename Imp::Traits::LocalBasisType::Traits> Interface;

    static Imp* clone(const Imp& imp)
    {
      return LocalFiniteElementCloneFactoryHelper<Imp, std::is_base_of<Interface, Imp>::value>::clone(imp);
    }

    static Imp* create()
    {
      return new Imp;
    }
  };



  // -----------------------------------------------------------------
  // Basis
  // -----------------------------------------------------------------

  /**
   * @brief class for wrapping a basis using the virtual interface
   *
   * The differentiation order of the traits T might be less than
   * the one in the traits of the implementation.
   *
   * @tparam T The LocalBasisTraits class
   * @tparam Imp LocalBasisInterface implementation
   */
  template<class T , class Imp>
  class LocalBasisVirtualImp
    : public virtual LocalBasisVirtualInterface<T>,
      public LocalBasisVirtualImp<typename LowerOrderLocalBasisTraits<T>::Traits,Imp>
  {
    template<class FEImp>
    friend class LocalFiniteElementVirtualImp;

    typedef LocalBasisVirtualImp<typename LowerOrderLocalBasisTraits<T>::Traits,Imp> Base;

  protected:

    //! constructor taking an implementation of the interface
    LocalBasisVirtualImp( const Imp &imp )
      : Base(imp)
    {}

  public:
    typedef T Traits;

    using Base::size;
    using Base::order;
    using Base::evaluateFunction;
    using Base::evaluateJacobian;
    using Base::evaluate;
    using Base::partial;


    //! @copydoc LocalBasisVirtualInterface::evaluate
    inline void evaluate(
      const std::array<int,Traits::diffOrder>& directions,
      const typename Traits::DomainType& in,
      std::vector<typename Traits::RangeType>& out) const
    {
      // Even for double virtualization it is save to call the template method
      // since the interface provides it redirecting to the virtual method
      // of the derived class
      //
      // Unfortunately not all compilers can determine Traits::diffOrder from
      // the type of the argument directions
      DUNE_NO_DEPRECATED_BEGIN
      impl_.template evaluate<Traits::diffOrder>(directions, in, out);
      DUNE_NO_DEPRECATED_END
    }

  protected:
    using Base::impl_;

  };


  /**
   * @brief class for wrapping a basis using the virtual interface
   *
   * This is the base class of all wrappers. It has differentiation order 0.
   *
   * @tparam Imp LocalBasisInterface implementation
   */
  template<class DF, int n, class D, class RF, int m, class R, class J, class Imp>
  class LocalBasisVirtualImp<LocalBasisTraits<DF,n,D,RF,m,R,J,0>, Imp>
    : public virtual LocalBasisVirtualInterface<LocalBasisTraits<DF,n,D,RF,m,R,J,0> >
  {
    template<class FEImp>
    friend class LocalFiniteElementVirtualImp;

  protected:

    //! constructor taking an implementation of the interface
    LocalBasisVirtualImp( const Imp &imp )
      : impl_(imp)
    {}

  public:
    typedef LocalBasisTraits<DF,n,D,RF,m,R,J,0> Traits;

    //! @copydoc LocalBasisVirtualInterface::size
    unsigned int size () const
    {
      return impl_.size();
    }

    //! @copydoc LocalBasisVirtualInterface::order
    unsigned int order () const
    {
      return impl_.order();
    }

    //! @copydoc LocalBasisVirtualInterface::evaluateFunction
    inline void evaluateFunction (const typename Traits::DomainType& in,
                                  std::vector<typename Traits::RangeType>& out) const
    {
      impl_.evaluateFunction(in,out);
    }

    //! @copydoc LocalBasisVirtualInterface::evaluateJacobian
    inline void evaluateJacobian(
      const typename Traits::DomainType& in,
      std::vector<typename Traits::JacobianType>& out) const
    {
      impl_.evaluateJacobian(in,out);
    }

    /** \brief Evaluate partial derivatives of any order of all shape functions
     * \param order Order of the partial derivatives, in the classic multi-index notation
     * \param in Position where to evaluate the derivatives
     * \param[out] out Return value: the desired partial derivatives
     */
    void partial(const std::array<unsigned int,Traits::dimDomain>& order,
                 const typename Traits::DomainType& in,
                 std::vector<typename Traits::RangeType>& out) const
    {
      impl_.partial(order,in,out);
    }

    //! @copydoc LocalBasisVirtualInterface::evaluate
    inline void evaluate(
      const std::array<int,Traits::diffOrder>& directions,
      const typename Traits::DomainType& in,
      std::vector<typename Traits::RangeType>& out) const
    {
      //      impl_.template evaluate<Traits::diffOrder> (directions, in,out);
      //      impl_.template evaluate<0> (directions, in,out);
      impl_.evaluateFunction(in,out);
    }

  protected:
    const Imp& impl_;
  };



  // -----------------------------------------------------------------
  // Interpolation
  // -----------------------------------------------------------------

  /**
   * @brief class for wrapping a local interpolation
   *        using the virtual interface
   *
   * @tparam DomainType domain type of the Dune::VirtualFunction to interpolate
   * @tparam RangeType range type of the Dune::VirtualFunction to interpolate
   * \tparam Imp LocalInterpolationVirtualInterface implementation
   */
  template<class DomainType, class RangeType, class Imp>
  class LocalInterpolationVirtualImp
    : public LocalInterpolationVirtualInterface< DomainType, RangeType >
  {
    template<class FEImp>
    friend class LocalFiniteElementVirtualImp;

    typedef LocalInterpolationVirtualInterface< DomainType, RangeType > Base;

  protected:

    //! constructor taking an implementation of the Dune::LocalInterpolationVirtualInterface
    LocalInterpolationVirtualImp( const Imp &imp)
      : impl_(imp) {}

  public:

    typedef typename Base::FunctionType FunctionType;

    typedef typename Base::CoefficientType CoefficientType;

    //! \copydoc LocalInterpolationVirtualInterface::interpolate
    virtual void interpolate (const FunctionType& f, std::vector<CoefficientType>& out) const
    {
      impl_.interpolate(f,out);
    }

  protected:
    const Imp& impl_;

  };



  // -----------------------------------------------------------------
  // Coefficients
  // -----------------------------------------------------------------

  /**
   * @brief class for wrapping local coefficients
   *        using the virtual interface
   *
   * @tparam Imp LocalCoefficientsInterface implementation
   */
  template<class Imp>
  class LocalCoefficientsVirtualImp
    : public LocalCoefficientsVirtualInterface
  {
    template<class FEImp>
    friend class LocalFiniteElementVirtualImp;

  protected:

    //! constructor taking an implementation of the Dune::LocalCoefficientsVirtualInterface
    LocalCoefficientsVirtualImp( const Imp &imp )
      : impl_(imp)
    {}

  public:

    //! @copydoc LocalCoefficientsVirtualInterface::size
    std::size_t size () const
    {
      return impl_.size();
    }

    //! @copydoc LocalCoefficientsVirtualInterface::localKey
    const LocalKey& localKey (std::size_t i) const
    {
      return impl_.localKey(i);
    }

  protected:
    const Imp& impl_;

  };



  // -----------------------------------------------------------------
  // Finite Element
  // -----------------------------------------------------------------

  /**
   * @brief class for wrapping a finite element using the virtual interface
   *
   * This automatically inherits the differentiation order of the wrapped
   * finite element and implements the corresponding interface
   *
   * @tparam Imp LocalBasisInterface implementation
   */
  template<class Imp>
  class LocalFiniteElementVirtualImp
    : public virtual LocalFiniteElementVirtualInterface<typename Imp::Traits::LocalBasisType::Traits>
  {
    typedef typename Imp::Traits::LocalBasisType::Traits T;
    typedef LocalFiniteElementVirtualInterface<T> Interface;

  public:
    typedef typename Interface::Traits Traits;

    //! @copydoc constructor taking a LocalFiniteElementVirtualInterface implementation
    LocalFiniteElementVirtualImp( const Imp &imp )
      : impl_(LocalFiniteElementCloneFactory<Imp>::clone(imp)),
        localBasisImp_(impl_->localBasis()),
        localCoefficientsImp_(impl_->localCoefficients()),
        localInterpolationImp_(impl_->localInterpolation())
    {}

    //! Default constructor.  Assumes that the implementation class is default constructible as well.
    LocalFiniteElementVirtualImp()
      : impl_(LocalFiniteElementCloneFactory<Imp>::create()),
        localBasisImp_(impl_->localBasis()),
        localCoefficientsImp_(impl_->localCoefficients()),
        localInterpolationImp_(impl_->localInterpolation())
    {}

    //! Copy contructor needed for deep copy
    LocalFiniteElementVirtualImp(const LocalFiniteElementVirtualImp& other)
      : impl_(LocalFiniteElementCloneFactory<Imp>::clone(*other.impl_)),
        localBasisImp_(impl_->localBasis()),
        localCoefficientsImp_(impl_->localCoefficients()),
        localInterpolationImp_(impl_->localInterpolation())
    {}

    ~LocalFiniteElementVirtualImp()
    {
      delete impl_;
    }

    //! \copydoc LocalFiniteElementVirtualInterface::localBasis
    const typename Traits::LocalBasisType& localBasis () const
    {
      return localBasisImp_;
    }

    //! \copydoc LocalFiniteElementVirtualInterface::localCoefficients
    const typename Traits::LocalCoefficientsType& localCoefficients () const
    {
      return localCoefficientsImp_;
    }

    //! \copydoc LocalFiniteElementVirtualInterface::localInterpolation
    const typename Traits::LocalInterpolationType& localInterpolation () const
    {
      return localInterpolationImp_;
    }

    /** \brief Number of shape functions in this finite element */
    unsigned int size () const
    {
      return impl_->size();
    }

    //! \copydoc LocalFiniteElementVirtualInterface::type
    const GeometryType type () const
    {
      return impl_->type();
    }

    /** \brief clone this wrapper
     *
     * This 'virtual copy constructor' is needed if you want to copy
     * the wrapper through the virtual interface.
     */
    virtual LocalFiniteElementVirtualImp<Imp>* clone() const
    {
      return new LocalFiniteElementVirtualImp<Imp>(*this);
    }

  protected:
    const Imp* impl_;

    /** \todo This needs to automatically change to C0LocalBasisBla... to work with C0 shape functions */
    const LocalBasisVirtualImp<T, typename Imp::Traits::LocalBasisType> localBasisImp_;
    const LocalCoefficientsVirtualImp<typename Imp::Traits::LocalCoefficientsType> localCoefficientsImp_;
    const LocalInterpolationVirtualImp<typename T::DomainType,
        typename T::RangeType,
        typename Imp::Traits::LocalInterpolationType> localInterpolationImp_;
  };
}
#endif
