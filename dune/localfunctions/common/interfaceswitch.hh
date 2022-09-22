// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

#ifndef DUNE_LOCALFUNCTIONS_COMMON_INTERFACESWITCH_HH
#define DUNE_LOCALFUNCTIONS_COMMON_INTERFACESWITCH_HH

#include <cstddef>
#include <memory>
#include <vector>

#include <dune/common/fmatrix.hh>
#include <dune/common/typetraits.hh>
#include <dune/common/shared_ptr.hh>

namespace Dune {

  //! \brief Switch for uniform treatment of finite element with either the
  //!        local or the global interface
  /**
   * \tparam FiniteElement Type of the finite element to handle.
   * \tparam Dummy         Dummy parameter for enable_if.  This must be left
   *                       at the default value of \c void.
   *
   * \note The local interface is detected by the presence of the type
   *       FiniteElement::Traits::LocalBasisType.
   */
  template<class FiniteElement, class Dummy = void>
  struct FiniteElementInterfaceSwitch {
    //! export the type of the basis
    typedef typename FiniteElement::Traits::Basis Basis;
    //! export the type of the interpolation
    typedef typename FiniteElement::Traits::Interpolation Interpolation;
    //! export the type of the coefficients
    typedef typename FiniteElement::Traits::Coefficients Coefficients;

    //! access basis
    static const Basis &basis(const FiniteElement& fe)
    { return fe.basis(); }
    //! access interpolation
    static const Interpolation &interpolation(const FiniteElement& fe)
    { return fe.interpolation(); }
    //! access coefficients
    static const Coefficients &coefficients(const FiniteElement& fe)
    { return fe.coefficients(); }

    //! Type for storing finite elements
    /**
     * Some algorithms use one variable to store (as a shared pointer)
     * a finite element and update that pointer while iterating
     * through the grid. This works well as long as there is only a
     * moderate number of different finite elements, which can be
     * stored somewhere and don't change for the duration of the
     * algorithm.  This is the case for most local finite elements,
     * since they exists in a finite number of variants.
     *
     * If the number of possible finite element realizations grows to
     * big, e.g. for global finite elements or also for p-adaptive
     * local finite elements, these are only created on the
     * fly. Therefore we need to store a copy.  Since finite elements
     * in general are not assignable, we either copy-construct or
     * move-construt them for each grid element we visit.
     *
     * To accommodate both interfaces, we store in a `shared_ptr`.
     * Different ways to initialize a possible, from an l-value
     * reference, an r-value reference or from a shared_ptr.
     *
     * For backwards compatibility we assume that an l-value reference
     * to a local finite element is persistent and that we can simply
     * store the pointer using `stackobject_to_shared_ptr`, while in
     * the case of global finite elements we always need to copy
     * construct using `make_shared`. If a local finite element is not
     * persistent, it should be passed in as an r-value
     * reference. Access to the finite element is done by simply
     * dereferencing the store in both cases.
     */
    typedef std::shared_ptr<const FiniteElement> Store;
    //! Store a finite element in the store.
    /**
     * For local finite elements this means storing the address of the passed
     * reference, for global finite element this means creating a new object
     * with allocation and copy-construction and storing that.
     */
    static void setStore(Store& store, const FiniteElement& fe)
    { store = std::make_shared<const FiniteElement>(fe); }
    //! Store a finite element in the store.
    static void setStore(Store& store, FiniteElement&& fe)
    { store = std::make_shared<const FiniteElement>(std::move(fe)); }
    //! Store a finite element in the store.
    static void setStore(Store& store, const Store& fe)
    { store = fe; }
  };

#ifndef DOXYGEN
  //! \brief Switch for uniform treatment of finite element with either the
  //!        local or the global interface
  template<class FiniteElement>
  struct FiniteElementInterfaceSwitch<
      FiniteElement,
      typename std::enable_if<AlwaysTrue<typename FiniteElement::Traits::
              LocalBasisType>::value>::type
      >
  {
    //! export the type of the basis
    typedef typename FiniteElement::Traits::LocalBasisType Basis;
    //! export the type of the interpolation
    typedef typename FiniteElement::Traits::LocalInterpolationType
    Interpolation;
    //! export the type of the coefficients
    typedef typename FiniteElement::Traits::LocalCoefficientsType Coefficients;

    //! access basis
    static const Basis &basis(const FiniteElement& fe)
    { return fe.localBasis(); }
    //! access interpolation
    static const Interpolation &interpolation(const FiniteElement& fe)
    { return fe.localInterpolation(); }
    //! access coefficients
    static const Coefficients &coefficients(const FiniteElement& fe)
    { return fe.localCoefficients(); }

    //! Type for storing finite elements
    typedef std::shared_ptr<const FiniteElement> Store;
    //! Store a finite element in the store.
    static void setStore(Store& store, const FiniteElement& fe)
    { store = stackobject_to_shared_ptr<const FiniteElement>(fe); }
    //! Store a finite element in the store.
    static void setStore(Store& store, FiniteElement&& fe)
    { store = std::make_shared<const FiniteElement>(std::move(fe)); }
    //! Store a finite element in the store.
    static void setStore(Store& store, const Store& fe)
    { store = fe; }
  };
#endif // !DOXYGEN

  //! Switch for uniform treatment of local and global basis classes
  /**
   * \tparam Basis Type of the basis to handle.

   * \tparam Dummy Dummy parameter for enable_if.  This must be left at the
   *               default value of \c void.
   *
   * We don't provide any uniform access to the types and constants pertaining
   * to the global domain.  Providing this would require the Geometry as
   * template parameter as well, and the user code can build them itself if it
   * needs them with the help of the geometry.  The omitted types are \c
   * DomainGlobal and \c Jacobian, the omitted constant is \c dimDomainGlobal.
   *
   * \note The local interface is assumed if the constant
   *       Basis::Traits::dimDomain exists and has a value greater than 0.
   */
  template<class Basis, class Dummy = void>
  struct BasisInterfaceSwitch {
    //! export field types of the coordinates
    typedef typename Basis::Traits::DomainField DomainField;
    //! export dimension of local coordinates
    static const std::size_t dimDomainLocal = Basis::Traits::dimDomainLocal;
    //! export vector type of the local coordinates
    typedef typename Basis::Traits::DomainLocal DomainLocal;

    //! export field type of the values
    typedef typename Basis::Traits::RangeField RangeField;
    //! export dimension of the values
    static const std::size_t dimRange = Basis::Traits::dimRange;
    //! export vector type of the values
    typedef typename Basis::Traits::Range Range;

    //! Compute global gradient for scalar valued bases
    /**
     * \param basis    The basis to get the derivatives from.
     * \param geometry The geometry to use to transform the derivatives (for a
     *                 local basis, unused in the case of a global basis).
     * \param xl       The local coordinates where to evaluate the gradient.
     * \param grad     The result (will be resized to the appropriate number
     *                 of entries.
     *
     * \note This make sense only for a scalar valued basis.
     */
    template<typename Geometry>
    static void gradient(const Basis& basis, const Geometry& geometry,
                         const DomainLocal& xl,
                         std::vector<FieldMatrix<RangeField, 1,
                                 Geometry::coorddimension> >& grad)
    {
      grad.resize(basis.size());
      basis.evaluateJacobian(xl, grad);
    }
  };

#ifndef DOXYGEN
  //! Switch for uniform treatment of local and global basis classes
  template<class Basis>
  struct BasisInterfaceSwitch<Basis,
                              typename std::enable_if<
                                AlwaysTrue<
                                  std::integral_constant<
                                    std::size_t,
                                    Basis::Traits::dimDomain
                                    >
                                  >::value
                                >::type
                              >
  {
    //! export field types of the coordinates
    typedef typename Basis::Traits::DomainFieldType DomainField;
    //! export dimension of local coordinates
    static const std::size_t dimDomainLocal = Basis::Traits::dimDomain;
    //! export vector type of the local coordinates
    typedef typename Basis::Traits::DomainType DomainLocal;

    //! export field type of the values
    typedef typename Basis::Traits::RangeFieldType RangeField;
    //! export dimension of the values
    static const std::size_t dimRange = Basis::Traits::dimRange;
    //! export vector type of the values
    typedef typename Basis::Traits::RangeType Range;

    //! Compute global gradient for scalar valued bases
    template<typename Geometry>
    static void gradient(const Basis& basis, const Geometry& geometry,
                         const DomainLocal& xl,
                         std::vector<FieldMatrix<RangeField, 1,
                                 Geometry::coorddimension> >& grad)
    {
      std::vector<typename Basis::Traits::JacobianType> lgrad(basis.size());
      basis.evaluateJacobian(xl, lgrad);

      const typename Geometry::JacobianInverseTransposed& jac =
        geometry.jacobianInverseTransposed(xl);

      grad.resize(basis.size());
      for(std::size_t i = 0; i < basis.size(); ++i)
        jac.mv(lgrad[i][0], grad[i][0]);
    }
  };
#endif // !DOXYGEN

} // namespace Dune

#endif // DUNE_LOCALFUNCTIONS_COMMON_INTERFACESWITCH_HH
