// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

#ifndef DUNE_LOCALFUNCTIONS_META_POWER_BASIS_HH
#define DUNE_LOCALFUNCTIONS_META_POWER_BASIS_HH

#include <numeric>
#include <cstddef>
#include <vector>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

namespace Dune {

  //! Meta-basis turning a scalar basis into vector-valued basis
  /**
   * @ingroup BasisImplementation
   *
   * \tparam Backend Type of basis to take the power of.
   * \tparam dimR    Power to raise the basis to.
   */
  template<class Backend, std::size_t dimR>
  class PowerBasis {
    static_assert(Backend::Traits::dimRange == 1,
                  "PowerBasis works only with scalar backends");

    // don't use a reference here so this class stays copyable
    const Backend *backend;

  public:
    //! types of domain and range
    struct Traits : public Backend::Traits
    {
      //! Dimension of the range values
      static const std::size_t dimRange = dimR;
      //! Type used for range values
      typedef FieldVector<typename Traits::RangeField, dimR> Range;

      //! Jacobian properties
      /**
       * \note The Jacobian should be some matrix type with \c dimRange x
       *       \c dimDomainGlobal components of type \c RangeField.
       */
      typedef FieldMatrix<typename Traits::RangeField, dimR,
          Traits::dimDomainGlobal> Jacobian;
    };

    //! Construct a PowerBasis
    /**
     * \param backend_ Backend basis object to construct this object from.
     *                 This object holds a reference to the backend object.
     *                 This reference is also copied when this object is
     *                 copied.
     */
    PowerBasis(const Backend &backend_) : backend(&backend_) { }

    //! Number of shape functions
    std::size_t size () const { return backend->size()*dimR; }
    //! Polynomial order of the shape functions for quadrature
    std::size_t order () const { return backend->order(); }

    //! Evaluate all shape functions at given position
    void evaluateFunction(const typename Traits::DomainLocal& in,
                          std::vector<typename Traits::Range>& out) const
    {
      std::vector<typename Backend::Traits::Range> backendValues;
      backend->evaluateFunction(in, backendValues);
      out.assign(size(), typename Traits::Range(0));
      for(std::size_t d = 0; d < dimR; ++d)
        for(std::size_t i = 0; i < backend->size(); ++i)
          out[d*backend->size()+i][d] = backendValues[i][0];
    }

    //! Evaluate Jacobian of all shape functions at given position
    void evaluateJacobian(const typename Traits::DomainLocal& in,
                          std::vector<typename Traits::Jacobian>& out) const
    {
      std::vector<typename Backend::Traits::Jacobian> backendValues;
      backend->evaluateJacobian(in, backendValues);
      out.assign(size(), typename Traits::Jacobian(0));
      for(std::size_t d = 0; d < dimR; ++d)
        for(std::size_t i = 0; i < backend->size(); ++i)
          out[d*backend->size()+i][d] = backendValues[i][0];
    }

    //! \brief Evaluate partial derivatives of all shape functions
    void partial (const std::array<unsigned int, Backend::Traits::dimDomainGlobal>& order,
                  const typename Traits::DomainLocal& in,         // position
                  std::vector<typename Traits::Range>& out) const      // return value
    {
      auto totalOrder = std::accumulate(order.begin(), order.end(), 0);
      if (totalOrder == 0) {
        evaluateFunction(in, out);
      } else {
        DUNE_THROW(NotImplemented, "Desired derivative order is not implemented");
      }
    }
  };

} // namespace Dune

#endif // DUNE_LOCALFUNCTIONS_META_POWER_BASIS_HH
