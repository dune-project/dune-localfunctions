// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_COMMON_LOCALINTERPOLATION_HH
#define DUNE_LOCALFUNCTIONS_COMMON_LOCALINTERPOLATION_HH

#include <functional>

#include <dune/common/concept.hh>



namespace Dune {

  namespace Impl {

    // Concept for function supporting f.evaluate(Domain, Range&)
    template<class Domain, class Range>
    struct FunctionWithEvaluate
    {
      template<class F>
      auto require(F&& f) -> decltype(
        f.evaluate(std::declval<Domain>(), std::declval<Range&>())
      );
    };

    // Concept for function supporting f(Domain)
    template<class Domain>
    struct FunctionWithCallOperator
    {
      template<class F>
      auto require(F&& f) -> decltype(
        f(std::declval<Domain>())
      );
    };

    // Create function supporting Range = f(Domain)
    // If the argument already does this, just forward it.
    template<class Domain, class F,
      std::enable_if_t<models<FunctionWithCallOperator<Domain>, F>(), int> = 0>
    decltype(auto) makeFunctionWithCallOperator(const F& f)
    {
      return f;
    }

    // Create function supporting Range = f(Domain)
    // If the argument does not support this, wrap it in a lambda
    template<class Domain, class F,
      std::enable_if_t<not models<FunctionWithCallOperator<std::decay_t<Domain> >, F>(), int> = 0>
#ifndef DUNE_DEPRECATED_INTERPOLATE_CHECK
    [[deprecated( "Passing functions only supporting 'f.evaluate(x,y)' to interpolate() is deprecated."
                  "Use functions supporting operator(), i.e. f(x) instead!")]]
#endif
    decltype(auto) makeFunctionWithCallOperator(const F& f)
    {
      return [&](auto&& x) {
        typename std::decay_t<F>::Traits::RangeType y;
        f.evaluate(x,y);
        return y;
      };
    }

  } // namespace Impl

} // namespace Dune
#endif
