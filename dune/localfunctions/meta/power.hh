// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

#ifndef DUNE_LOCALFUNCTIONS_META_POWER_HH
#define DUNE_LOCALFUNCTIONS_META_POWER_HH

#include <cstddef>
#include <memory>

#include <dune/geometry/type.hh>

#include <dune/localfunctions/meta/power/basis.hh>
#include <dune/localfunctions/meta/power/coefficients.hh>
#include <dune/localfunctions/meta/power/interpolation.hh>

namespace Dune {

  //! \brief Meta-finite element turning a scalar finite element into
  //!        vector-valued one
  /**
   * @ingroup FiniteElementImplementation
   *
   * \tparam Backend Type of finite element to take the power of.
   * \tparam dimR    Power to raise the finite element to.
   */
  template<class Backend, std::size_t dimR>
  class PowerFiniteElement {
  public:
    //! types of component objects
    struct Traits {
      //! type of the Basis
      typedef PowerBasis<typename Backend::Traits::Basis, dimR> Basis;
      //! type of the Coefficients
      typedef PowerCoefficients Coefficients;
      //! type of the Interpolation
      typedef PowerInterpolation<typename Backend::Traits::Interpolation,
          typename Basis::Traits> Interpolation;
    };
  private:
    std::shared_ptr<const Backend> backend;
    typename Traits::Basis basis_;
    typename Traits::Coefficients coefficients_;
    typename Traits::Interpolation interpolation_;

  public:
    //! Construct a finite element
    /**
     * \note With this constructor a copy of the backend finite element is
     *       stored in this object.
     */
    PowerFiniteElement(const Backend &backend_) :
      backend(new Backend(backend_)),
      basis_(backend->basis()),
      coefficients_(backend->coefficients(), dimR),
      interpolation_(backend->interpolation())
    { }

    //! Construct a finite element
    /**
     * \note With this constructor ownership of the backend finite element is
     *       determined by the shared_ptr.
     */
    PowerFiniteElement(const std::shared_ptr<const Backend> &backendSPtr) :
      backend(backendSPtr),
      basis_(backend->basis()),
      coefficients_(backend->coefficients(), dimR),
      interpolation_(backend->interpolation())
    { }

    //! Extract basis of this finite element
    /**
     * The returned lvalue must have a lifetime at least as long as the finite
     * element object it was acquired from.
     */
    const typename Traits::Basis& basis() const { return basis_; }
    //! Extract coefficients of this finite element
    /**
     * The returned lvalue must have a lifetime at least as long as the finite
     * element object it was acquired from.
     */
    const typename Traits::Coefficients& coefficients() const
    { return coefficients_; }
    //! Extract interpolation of this finite element
    /**
     * The returned lvalue must have a lifetime at least as long as the finite
     * element object it was acquired from.
     */
    const typename Traits::Interpolation& interpolation() const
    { return interpolation_; }
    //! Extract geometry type of this finite element
    GeometryType type() const { return backend->type(); }
  };

  //! \brief Factory for meta-finite elements turning scalar finite elements
  //!        into vector-valued ones
  /**
   * \implements FiniteElementFactory
   * \ingroup FiniteElementFactoryImplementation
   *
   * \tparam BackendFiniteElement Type of finite element to take the power of.
   * \tparam dimR                 Power to raise the finite element to.
   */
  template<class BackendFiniteElement, std::size_t dimR>
  class PowerFiniteElementFactory
  {
  public:
    //! Type of the finite element
    typedef PowerFiniteElement<BackendFiniteElement, dimR> FiniteElement;

    //! create a finite element
    /**
     * \note With this overload of make() the backend finite element is copied
     *       into the created object.
     */
    const FiniteElement make(const BackendFiniteElement &backend) const
    { return FiniteElement(backend); }
    //! create a finite element
    /**
     * \note With this overload of make() ownership of the backend finite
     *       element is determined by the shared_ptr.
     */
    const FiniteElement
    make(const std::shared_ptr<const BackendFiniteElement> &backendSPtr) const
    { return FiniteElement(backendSPtr); }

  };

} // namespace Dune

#endif // DUNE_LOCALFUNCTIONS_META_POWER_HH
