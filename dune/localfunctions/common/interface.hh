// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

#ifndef DUNE_LOCALFUNCTIONS_INTERFACE_HH
#define DUNE_LOCALFUNCTIONS_INTERFACE_HH

#ifndef HEADERCHECK
#error This header exists for documentation purposes only and should never be included directly.
#endif

#include <array>
#include <cstddef>
#include <vector>


#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/localkey.hh>

namespace Dune {

  //! Interface for global-valued finite elements
  class FiniteElementInterface
  {
    struct ImplementationDefined;

  public:
    //! types of component objects
    /**
     * \note This may be a typedef instead of a member class.
     */
    struct Traits
    {
      //! type of the Basis
      /**
       * Should be an implementation of BasisInterface
       *
       * \note May be an inline class instead of a typedef.
       */
      typedef ImplementationDefined Basis;
      //! type of the Coefficients
      /**
       * Should be an implementation of CoefficientsInterface
       *
       * \note May be an inline class instead of a typedef.
       */
      typedef ImplementationDefined Coefficients;
      //! type of the Interpolation
      /**
       * Should be an implementation of InterpolationInterface
       *
       * \note May be an inline class instead of a typedef.
       */
      typedef ImplementationDefined Interpolation;
    };

    //! Construct a finite element
    /**
     * \note The arguments of the constructor are implementation specific.  In
     *       fact, finite element implementations are not required to be
     *       constructible by the user at all (except for copy-construction).
     *       The official way to construct a finite element is to use its
     *       factory.
     */
    FiniteElementInterface(...);
    //! Finite elements are CopyConstructible
    FiniteElementInterface(const FiniteElementInterface&);

    //! Extract basis of this finite element
    /**
     * The returned lvalue must have a lifetime at least as long as the finite
     * element object it was acquired from.
     */
    const Traits::Basis& basis() const;
    //! Extract coefficients of this finite element
    /**
     * The returned lvalue must have a lifetime at least as long as the finite
     * element object it was acquired from.
     */
    const Traits::Coefficients& coefficients() const;
    //! Extract interpolation of this finite element
    /**
     * The returned lvalue must have a lifetime at least as long as the finite
     * element object it was acquired from.
     */
    const Traits::Interpolation& interpolation() const;
    //! Extract geometry type of this finite element
    GeometryType type() const;
  };

  //! Factory interface for global-valued finite elements
  /**
   * The main purpose of the factory class is to provide a concept for
   * caching.  Take for instance a global-valued finite element that wraps a
   * local finite element.  The local finite element will typically have very
   * few variants, and the global-valued finite element will just apply a
   * geometric transformation to the derivatives.  The wrapped local finite
   * elements can be stored inside the factory and any global-valued finite
   * elements created by the factory just contain references or pointers.
   * This way the local finite elements don't need to be created anew for each
   * global-valued finite element.
   *
   * The other purpose is to semi-standardize the interface used to actually
   * create finite elements.  "Semi" because the information needed to create
   * an actual global-valued finite element will vary between finite element
   * types.  There are however certain types of information that are needed by
   * a larger subset of all available finite elements, so it makes sense to
   * define a common encoding for these types of information.  On the other
   * hand this information is often expensive to obtain, so it makes sense to
   * only provide it when it is actually needed.
   */
  template<class Geometry, class VertexOrder>
  class FiniteElementFactoryInterface
  {
    struct ImplementationDefined;

  public:
    //! Type of the finite element
    /**
     * Should be an implementation of FiniteElementInterface
     *
     * \note May be an inline class instead of a typedef.
     */
    typedef ImplementationDefined FiniteElement;

    //! Construct a finite element factory
    /**
     * \note The arguments of the constructor are implementation specific.
     */
    FiniteElementFactoryInterface(...);

    /** \name Finite element creation methods
     *
     * Each finite element factory implementation must provide at least one of
     * these methods.  The signatures may be extended by additional
     * parameters, but the parameters that are specified here should be given
     * first and in the order specified here.
     *
     * The return value of these functions is suitable for binding to a const
     * reference -- it will either be an rvalue (in which case binding to a
     * const reference will create a copy whose lifetime is the same as the
     * reference itself), or it will be an lvalue (in which case the factory
     * must guarantee that it will be valid until the factory is destroyed or
     * something else happens that explicitly invalidates all created finite
     * elements).
     *
     * In any case, since global-valued finite element objects are
     * copy-constructible, it is also possible to use the returned value to
     * initialize a finite element object instead of a const reference.
     */
    //! \{

    //! create a finite element from a geometry and a vertex ordering
    const FiniteElement make(const Geometry&, const VertexOrder&, ...);
    //! create a finite element from a geometry
    const FiniteElement make(const Geometry&, ...);
    //! create a finite element from a vertex ordering
    const FiniteElement make(const VertexOrder&, ...);
    //! create a finite element from a geometry type
    /**
     * \note This signature should only be used when only the geometry type
     *       but not the full geometry or vertex ordering are needed.
     */
    const FiniteElement make(const GeometryType&, ...);
    //! create a finite element
    const FiniteElement make(...);

    //! \}

  };

  //! Interface for global-valued shape functions
  class BasisInterface
  {
    struct ImplementationDefined;
    constexpr static int implementationDefined = 42;

  public:
    //! types of domain and range
    /**
     * \nosubgrouping
     *
     * \note This may be a typedef instead of a member class.
     */
    struct Traits
    {
      //! \name Domain properties (local and global)
      //! \{

      //! Field type of the domain
      typedef ImplementationDefined DomainFieldType;

      //! \brief dimension of the domain
      constexpr static int dimDomain = implementationDefined;

      //! Type used for coordinate vectors in the domain
      typedef ImplementationDefined DomainType;

      //! \}

      //! \name Range properties (global range only)
      //! \{

      //! Field type of the range
      typedef ImplementationDefined RangeFieldType;

      //! \brief dimension of the range
      constexpr static int dimRange = implementationDefined;

      //! Type used for range values
      typedef ImplementationDefined RangeType;

      //! \}

      //! Jacobian properties
      /**
       * \note The Jacobian should be some matrix type with \c dimRange x
       *       \c dimDomain components of type \c RangeFieldType.
       */
      typedef ImplementationDefined Jacobian;
    };

    //! Number of shape functions
    std::size_t size () const;
    //! Polynomial order of the shape functions for quadrature
    std::size_t order () const;

    //! Evaluate all shape functions at given position
    void evaluateFunction(const Traits::DomainType& in,
                          std::vector<Traits::RangeType>& out) const;

    //! Evaluate Jacobian of all shape functions at given position
    void evaluateJacobian(const Traits::DomainType& in,
                          std::vector<Traits::Jacobian>& out) const;

    /** \brief Evaluate partial derivatives of any order of all shape functions
     * \param order Order of the partial derivatives, in the classic multi-index notation
     * \param in Position where to evaluate the derivatives
     * \param[out] out Return value: the desired partial derivatives
     */
    void partial(const std::array<unsigned int,Traits::dimDomain>& order,
                 const typename Traits::DomainType& in,
      std::vector<typename Traits::RangeType>& out) const;
  };

  //! Interface for global-valued interpolation
  struct InterpolationInterface
  {
    //! Export basis traits
    /**
     * This should be the traits class of the corresponding basis.
     */
    typedef BasisInterface::Traits Traits;

    //! Determine coefficients interpolating a given function
    /**
     * \param f   An object supporting the expression \c f.evaluate(x,y),
     *            where \c x is of type \c Traits::DomainLocal and \c y of the
     *            type \c Traits::Range.  When \c f.evaluate(x,y) is
     *            evaluated, \c x will be a local coordinate , and the
     *            expression should set \c y to the function value at that
     *            position.  The initial value of \c y should not be used.
     * \param out Vector where to store the interpolated coefficients.
     */
    template<typename F, typename C>
    void interpolate (const F& f, std::vector<C>& out) const;
  };

  //! Interface for global-valued coefficients
  /**
   * \note This interface is listed separately only to keep it together with
   *       the other global-valued interfaces.  It is identical to the
   *       interface for local coefficients.
   */
  struct CoefficientsInterface
  {
    //! number of coefficients
    std::size_t size() const;

    //! get i'th index
    const LocalKey& localKey(std::size_t i) const;
  };
}
#endif // DUNE_LOCALFUNCTIONS_INTERFACE_HH
