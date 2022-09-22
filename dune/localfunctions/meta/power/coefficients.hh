// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

#ifndef DUNE_LOCALFUNCTIONS_META_POWER_COEFFICIENTS_HH
#define DUNE_LOCALFUNCTIONS_META_POWER_COEFFICIENTS_HH

#include <cstddef>
#include <vector>

#include <dune/localfunctions/common/localkey.hh>

namespace Dune {

  //! \brief Meta-coefficients turning a scalar coefficients into
  //!        vector-valued coefficients
  /**
   * \nosubgrouping
   * \implements CoefficientsInterface
   */
  class PowerCoefficients {
    std::vector<LocalKey> keys;

  public:
    //! Construct a PowerCoefficients object
    /**
     * \param backend The backend coeficients object to raise to a power.
     * \param power   Power to raise the backend object to.
     *
     * The LocalKeys of the backend coefficients are copied into internal
     * storage.  The index member of each LocalKey is modified to keep them
     * unique for instances of different power.
     */
    template<class Backend>
    PowerCoefficients(const Backend &backend, std::size_t power) :
      keys(backend.size()*power)
    {
      for(std::size_t i = 0; i < backend.size(); ++i) {
        const LocalKey &k = backend.localKey(i);
        for(std::size_t d = 0; d < power; ++d)
          keys[i+d*backend.size()] =
            LocalKey(k.subEntity(), k.codim(), power*k.index() + d);
      }
    }
    //! number of coefficients
    inline std::size_t size() const { return keys.size(); }

    //! get i'th index
    inline const LocalKey& localKey(std::size_t i) const { return keys[i]; }
  };

} // namespace Dune

#endif // DUNE_LOCALFUNCTIONS_META_POWER_COEFFICIENTS_HH
