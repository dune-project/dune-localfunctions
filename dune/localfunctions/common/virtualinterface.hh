// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_COMMON_VIRTUALINTERFACE_HH
#define DUNE_LOCALFUNCTIONS_COMMON_VIRTUALINTERFACE_HH

#include <type_traits>
#include <array>
#include <vector>
#include <functional>

#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/localbasis.hh>
#include <dune/localfunctions/common/localinterpolation.hh>
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
   * @brief Return a proper base class for functions to use with LocalInterpolation.
   *
   * @tparam FE A FiniteElement type
   *
   * \deprecated
   * This class is deprecated.
   * To keep this traits class working it exports a simple
   * look-a-like of the old Dune::Function base class.
   * However, you should stop using this and pass functions with
   * plain operator() interface to interpolate() from now on.
   */
  template<class FE>
  class
  [[deprecated("Dune::LocalFiniteElementFunctionBase is deprecated after Dune 2.7. You can now pass functions providing operator() to interpolate.")]]
  LocalFiniteElementFunctionBase
  {
    typedef typename FE::Traits::LocalBasisType::Traits::DomainType Domain;
    typedef typename FE::Traits::LocalBasisType::Traits::RangeType Range;

    // Hack: Keep a copy of Dune::Function here. This allows to avoid depending
    // on the deprecated dune-common header while still keeping the LocalFiniteElementFunctionBase
    // mechanism working during its deprecation period.
    class FunctionBaseDummy
    {
    public:

      using RangeType = Range;
      using DomainType = Domain;

      struct Traits
      {
        using RangeType = Range;
        using DomainType = Domain;
      };

      void evaluate(const DomainType& x, RangeType& y) const;
    };

  public:

    using VirtualFunctionBase = FunctionBaseDummy;
    using FunctionBase = FunctionBaseDummy;

    /** \brief Base class type for functions to use with LocalInterpolation
     *
     * This is just a dummy providing the old typedefs.
     * interface and Function base class otherwise.
     */
    using type = FunctionBaseDummy;
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

    //! type of function to interpolate
    using FunctionType = std::function<RangeType(DomainType)>;

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

    //! type of function to interpolate
    using FunctionType = std::function<RangeType(DomainType)>;

    //! type of the coefficient vector in the interpolate method
    typedef typename RangeType::field_type CoefficientType;


    virtual ~LocalInterpolationVirtualInterface() {}

    // This method is only noted again for to make the documentation complete.

    /** \brief determine coefficients interpolating a given function
     *
     * This is the pure virtual method taking a VirtualFunction.
     *
     * \param[in]  f   Function instance used to interpolate.
     * \param[out] out Resulting coefficients vector.
     */
    virtual void interpolate (const FunctionType& f, std::vector<CoefficientType>& out) const = 0;

    /** \brief determine coefficients interpolating a given function
     *
     * \param[in]  ff   Function instance used to interpolate.
     * \param[out] out Resulting coefficients vector.
     */
    template<class F,
      std::enable_if_t<not std::is_base_of<FunctionType, F>::value, int> = 0>
    void interpolate (const F& ff, std::vector<CoefficientType>& out) const
    {
      const auto& f = Impl::makeFunctionWithCallOperator<DomainType>(ff);

      const LocalInterpolationVirtualInterfaceBase<DomainType, RangeType>& asBase = *this;
      asBase.interpolate(FunctionType(std::cref(f)),out);
    }

    /** \brief determine coefficients interpolating a given function
     *
     * \param[in]  ff   Function instance used to interpolate.
     * \param[out] out Resulting coefficients vector.
     */
    template<class F, class C>
    void interpolate (const F& ff, std::vector<C>& out) const
    {
      const auto& f = Impl::makeFunctionWithCallOperator<DomainType>(ff);

      std::vector<CoefficientType> outDummy;
      const LocalInterpolationVirtualInterfaceBase<DomainType, RangeType>& asBase = *this;
      asBase.interpolate(FunctionType(std::cref(f)),outDummy);
      out.resize(outDummy.size());
      for(typename std::vector<CoefficientType>::size_type i=0; i<outDummy.size(); ++i)
        out[i] = outDummy[i];
    }
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


  /**
   * @brief virtual base class for local finite elements with functions
   *
   * This class defines the same interface using pure virtual methods.
   */
  template<class T>
  class LocalFiniteElementVirtualInterface
  {
    using LocalBasisTraits = T;
  public:
    typedef LocalFiniteElementTraits<
        LocalBasisVirtualInterface<LocalBasisTraits>,
        LocalCoefficientsVirtualInterface,
        LocalInterpolationVirtualInterface<
            typename LocalBasisTraits::DomainType,
            typename LocalBasisTraits::RangeType> > Traits;

    virtual ~LocalFiniteElementVirtualInterface() {}

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

    virtual LocalFiniteElementVirtualInterface<T>* clone() const = 0;
  };
}
#endif
