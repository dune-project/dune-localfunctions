// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LOCALFUNCTIONS_COMMON_VIRTUALINTERFACE_HH
#define DUNE_LOCALFUNCTIONS_COMMON_VIRTUALINTERFACE_HH

#include <array>

#include <dune/common/function.hh>

#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/localbasis.hh>
#include <dune/localfunctions/common/localkey.hh>
#include <dune/localfunctions/common/localfiniteelementtraits.hh>

namespace Dune
{

  // forward declaration needed by the helper traits
  template<class DomainType, class RangeType>
  class LocalInterpolationVirtualInterface;

  // -----------------------------------------------------------------
  // Helper traits classes
  // -----------------------------------------------------------------

  /**
   * @brief Construct LocalBasisTraits with fixed diff order
   *
   * @tparam T A LocalBasisTraits class
   * @tparam order The differentiation order
   */
  template<class T, int order>
  struct FixedOrderLocalBasisTraits
  {
    //! The LocalBasisTraits specified order
    typedef LocalBasisTraits<
        typename T::DomainFieldType,
        T::dimDomain,
        typename T::DomainType,
        typename T::RangeFieldType,
        T::dimRange,
        typename T::RangeType,
        typename T::JacobianType,
        order> Traits;
  };

  /**
   * @brief Return a proper base class for functions to use with LocalInterpolation.
   *
   * @tparam FE A FiniteElement type
   */
  template<class FE>
  class LocalFiniteElementFunctionBase
  {
    typedef typename FE::Traits::LocalBasisType::Traits::DomainType DomainType;
    typedef typename FE::Traits::LocalBasisType::Traits::RangeType RangeType;

    typedef LocalInterpolationVirtualInterface<DomainType, RangeType> Interface;
    typedef typename FE::Traits::LocalInterpolationType Implementation;

  public:

    typedef VirtualFunction<DomainType, RangeType> VirtualFunctionBase;
    typedef Function<const DomainType&, RangeType&> FunctionBase;

    /** \brief Base class type for functions to use with LocalInterpolation
     *
     * This is the VirtualFunction interface class if FE implements the virtual
     * interface and Function base class otherwise.
     */
    typedef typename std::conditional<std::is_base_of<Interface, Implementation>::value, VirtualFunctionBase, FunctionBase>::type type;
  };



  // -----------------------------------------------------------------
  // Basis
  // -----------------------------------------------------------------

  /**
   * @brief virtual base class for a local basis
   *
   * Provides the local basis interface with pure virtual methods.
   * This class defines the interface using pure virtual methods.
   */
  template<class T>
  class LocalBasisVirtualInterface
  {
  public:
    using Traits = T;


    virtual ~LocalBasisVirtualInterface() {}

    //! \brief Number of shape functions
    virtual unsigned int size () const = 0;

    //! \brief Polynomial order of the shape functions
    virtual unsigned int order () const = 0;

    /** \brief Evaluate all basis function at given position
     *
     * Evaluates all shape functions at the given position and returns
     * these values in a vector.
     */
    virtual void evaluateFunction (const typename Traits::DomainType& in,
                                   std::vector<typename Traits::RangeType>& out) const = 0;

    /** \brief Evaluate jacobian of all shape functions at given position.
     *
     * out[i][j][k] is \f$\partial_k \hat\phi_j^i \f$, when \f$\hat\phi^i \f$ is the
     *		i'th shape function.
     *
     * \param [in]  in  The position where evaluated
     * \param [out] out The result
     */
    virtual void evaluateJacobian(const typename Traits::DomainType& in,         // position
                                  std::vector<typename Traits::JacobianType>& out) const = 0;

    /** \brief Evaluate partial derivatives of any order of all shape functions
     * \param order Order of the partial derivatives, in the classic multi-index notation
     * \param in Position where to evaluate the derivatives
     * \param[out] out Return value: the desired partial derivatives
     */
    virtual void partial(const std::array<unsigned int,Traits::dimDomain>& order,
                         const typename Traits::DomainType& in,
                         std::vector<typename Traits::RangeType>& out) const = 0;
  };



  // -----------------------------------------------------------------
  // Interpolation
  // -----------------------------------------------------------------

  /**
   * @brief virtual base class for a local interpolation
   *
   * This class defines the interface using pure virtual methods.
   * In applications you should use the derived class
   * LocalInterpolationVirtualInterface that also
   * contains a interpolate method where the function type
   * is a template parameter.
   *
   * This template method cannot be defined in the same
   * class as the virtual method. Otherwise name resolution fails.
   */
  template<class DomainType, class RangeType>
  class LocalInterpolationVirtualInterfaceBase
  {
  public:

    //! type of virtual function to interpolate
    typedef Dune::VirtualFunction<DomainType, RangeType> FunctionType;

    //! type of the coefficient vector in the interpolate method
    typedef typename RangeType::field_type CoefficientType;

    virtual ~LocalInterpolationVirtualInterfaceBase() {}

    /** \brief determine coefficients interpolating a given function
     *
     * This is the pure virtual method taking a VirtualFunction.
     *
     * \param[in]  f   Function instance used to interpolate.
     * \param[out] out Resulting coefficients vector.
     */
    virtual void interpolate (const FunctionType& f, std::vector<CoefficientType>& out) const = 0;
  };

  /**
   * @brief virtual base class for a local interpolation
   *
   * This class defines the interface using pure virtual methods.
   * It also contains the interpolate method with function
   * type as template parameter.
   */
  template<class DomainType, class RangeType>
  class LocalInterpolationVirtualInterface
    : public LocalInterpolationVirtualInterfaceBase<DomainType, RangeType>
  {
  public:

    //! type of virtual function to interpolate
    typedef Dune::VirtualFunction<DomainType, RangeType> FunctionType;

    //! type of the coefficient vector in the interpolate method
    typedef typename RangeType::field_type CoefficientType;


    virtual ~LocalInterpolationVirtualInterface() {}

    // This method is only notet again for to make the documentation complete.

    /** \brief determine coefficients interpolating a given function
     *
     * This is the pure virtual method taking a VirtualFunction.
     *
     * \param[in]  f   Function instance used to interpolate.
     * \param[out] out Resulting coefficients vector.
     */
    virtual void interpolate (const FunctionType& f, std::vector<CoefficientType>& out) const = 0;

    //! \copydoc LocalInterpolationVirtualInterfaceBase::interpolate
    //! This uses the pure virtual method by wrapping the template argument into a VirtualFunction
    template<class F>
    void interpolate (const F& f, std::vector<CoefficientType>& out) const
    {
      const LocalInterpolationVirtualInterfaceBase<DomainType, RangeType>& asBase = *this;
      asBase.interpolate(VirtualFunctionWrapper<F>(f),out);
    }

    template<class F, class C>
    void interpolate (const F& f, std::vector<C>& out) const
    {
      std::vector<CoefficientType> outDummy;
      const LocalInterpolationVirtualInterfaceBase<DomainType, RangeType>& asBase = *this;
      asBase.interpolate(VirtualFunctionWrapper<F>(f),outDummy);
      out.resize(outDummy.size());
      for(typename std::vector<CoefficientType>::size_type i=0; i<outDummy.size(); ++i)
        out[i] = outDummy[i];
    }

  private:

    template <typename F>
    struct VirtualFunctionWrapper
      : public FunctionType
    {
    public:
      VirtualFunctionWrapper(const F &f)
        : f_(f)
      {}

      virtual ~VirtualFunctionWrapper() {}

      virtual void evaluate(const DomainType& x, RangeType& y) const
      {
        f_.evaluate(x,y);
      }

      const F &f_;
    };
  };



  // -----------------------------------------------------------------
  // Coefficients
  // -----------------------------------------------------------------

  /**
   * @brief virtual base class for local coefficients
   *
   * This class defines the interface using pure virtual methods.
   */
  class LocalCoefficientsVirtualInterface
  {
  public:

    virtual ~LocalCoefficientsVirtualInterface() {}

    //! number of coefficients
    virtual std::size_t size () const = 0;

    //! get i'th index
    const virtual LocalKey& localKey (std::size_t i) const = 0;

  };



  // -----------------------------------------------------------------
  // Finite Element
  // -----------------------------------------------------------------


  /*
   * This class is necessary to canonicalize the diffOrder
   * in the LocalDerivativeTraits to the value zero. Doing
   * this via the template alias LocalFiniteElementVirtualInterface
   * has the consequence that you get exactly the same interface
   * class regardless of the diffOrder in the passed LocalBasisTraits.
   */
  template<class T>
  class LocalFiniteElementVirtualInterfaceImp
  {
    using LocalBasisTraits = T;
  public:
    typedef LocalFiniteElementTraits<
        LocalBasisVirtualInterface<LocalBasisTraits>,
        LocalCoefficientsVirtualInterface,
        LocalInterpolationVirtualInterface<
            typename LocalBasisTraits::DomainType,
            typename LocalBasisTraits::RangeType> > Traits;

    virtual ~LocalFiniteElementVirtualInterfaceImp() {}

    //! \copydoc LocalFiniteElementVirtualInterface::localBasis
    virtual const typename Traits::LocalBasisType& localBasis () const = 0;

    //! \copydoc LocalFiniteElementVirtualInterface::localCoefficients
    virtual const typename Traits::LocalCoefficientsType& localCoefficients () const = 0;

    //! \copydoc LocalFiniteElementVirtualInterface::localInterpolation
    virtual const typename Traits::LocalInterpolationType& localInterpolation () const = 0;

    //! \copydoc LocalFiniteElementVirtualInterface::size
    virtual unsigned int size () const = 0;

    //! \copydoc LocalFiniteElementVirtualInterface::type
    virtual const GeometryType type () const = 0;

    virtual LocalFiniteElementVirtualInterfaceImp<T>* clone() const = 0;
  };

  /**
   * @brief virtual base class for local finite elements with functions
   *
   * This class defines the same interface using pure virtual methods.
   */
  template<class T>
  using LocalFiniteElementVirtualInterface = LocalFiniteElementVirtualInterfaceImp<typename FixedOrderLocalBasisTraits<T,0>::Traits>;
}
#endif
