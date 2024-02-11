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
    // This functions returns the passed in reference.
    template<class Domain, class F,
      std::enable_if_t<models<FunctionWithCallOperator<Domain>, F>(), int> = 0>
    [[deprecated( "The utility function makeFunctionWithCallOperator() is deprecated and will be removed after 2.10."
                  "Downstream modules no longer need to call this function since interpolate() no-longer supports non-callable functions.")]]
    decltype(auto) makeFunctionWithCallOperator(const F& f)
    {
      return f;
    }

  } // namespace Impl

} // namespace Dune
#endif
