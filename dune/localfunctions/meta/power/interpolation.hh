// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

#ifndef DUNE_LOCALFUNCTIONS_META_POWER_INTERPOLATION_HH
#define DUNE_LOCALFUNCTIONS_META_POWER_INTERPOLATION_HH

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <vector>
#include <dune/localfunctions/common/localinterpolation.hh>

namespace Dune {

  //! \brief Meta-interpolation turning a scalar interpolation into
  //!        vector-valued interpolation
  /**
   * \tparam Backend     Type of the scalar interpolation.
   * \tparam BasisTraits Traits type of the corresponding PowerBasis.
   *
   * \nosubgrouping
   */
  template<class Backend, class BasisTraits>
  class PowerInterpolation {
    static_assert(Backend::Traits::dimRange == 1,
                  "PowerInterpolation  works only with scalar backends");

    const Backend *backend;

  public:
    //! Export basis traits
    typedef BasisTraits Traits;

    //! Construct a PowerInterpolation
    /**
     * \param backend_ Backend interpolation object to construct this object
     *                 from.  This object holds a reference to the backend
     *                 object.  This reference is also copied when this object
     *                 is copied.
     */
    PowerInterpolation(const Backend &backend_) : backend(&backend_) { }

  private:
    template<class F>
    class ComponentEvaluator
    {
      const F &f;
      std::size_t comp;

    public:
      ComponentEvaluator(const F &f_, std::size_t comp_) :
        f(f_), comp(comp_)
      { }

      typename Backend::Traits::Range operator()(const typename Backend::Traits::DomainLocal &x) const
      {
        typename Traits::Range fy = f(x);
        typename Backend::Traits::Range y;
        y[0] = fy[comp];
        return y;
      }
    };

  public:
    //! Determine coefficients interpolating a given function
    /**
     * \param ff  An object supporting the expression \c ff.evaluate(x,y),
     *            where \c x is of type \c Traits::DomainLocal and \c y of the
     *            type \c Traits::Range.  When \c f.evaluate(x,y) is
     *            evaluated, \c x will be a local coordinate, and the
     *            expression should set \c y to the function value at that
     *            position.  The initial value of \c y should not be used.
     * \param out Vector where to store the interpolated coefficients.
     */
    template<typename F, typename C>
    void interpolate(const F& ff, std::vector<C>& out) const {

      auto&& f = Impl::makeFunctionWithCallOperator<typename Backend::Traits::DomainLocal>(ff);


      out.clear();
      std::vector<C> cout;
      for(std::size_t d = 0; d < Traits::dimRange; ++d) {
        // When dropping support for `evaluate()` we can simply use a lambda
        // instead of ComponentEvaluator. But changing this now would break
        // PowerInterpolation for FE-implementation outside of dune-localfunctions
        // which may not have been adjusted so far.
        backend->interpolate(ComponentEvaluator<std::decay_t<decltype(f)>>(f, d), cout);
        if(d == 0)
          out.resize(cout.size()*Traits::dimRange);
        // make sure the size of cout does not change surprisingly
        assert(out.size() == cout.size()*Traits::dimRange);
        std::copy(cout.begin(), cout.end(), out.begin() + d*cout.size());
      }
    }
  };

} // namespace Dune

#endif // DUNE_LOCALFUNCTIONS_META_POWER_INTERPOLATION_HH
