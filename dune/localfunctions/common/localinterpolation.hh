// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LOCALFUNCTIONS_COMMON_LOCALINTERPOLATION_HH
#define DUNE_LOCALFUNCTIONS_COMMON_LOCALINTERPOLATION_HH

#include <functional>

#include <dune/common/concept.hh>
#include <dune/common/function.hh>



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

    // Create function supporting f.evaluate(Domain, Range&)
    // If the argument already does this, just forward it.
    template<class Domain, class Range, class F,
      std::enable_if_t<models<FunctionWithEvaluate<Domain, Range>, F>(), int> = 0>
    decltype(auto) makeFunctionWithEvaluate(const F& f)
    {
      return f;
    }

    // Create function supporting f.evaluate(Domain, Range&)
    // If the argument does not support this, wrap it as VirtualFunction
    template<class Domain, class Range, class F,
      std::enable_if_t<not models<FunctionWithEvaluate<Domain, Range>, F>(), int> = 0>
    decltype(auto) makeFunctionWithEvaluate(const F& f)
    {
      return makeVirtualFunction<Domain, Range>(std::cref(f));
    }

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
      std::enable_if_t<not models<FunctionWithCallOperator<Domain>, F>(), int> = 0>
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
