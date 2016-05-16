// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LOCALFUNCTIONS_COMMON_STATICFORLOOP_HH
#define DUNE_LOCALFUNCTIONS_COMMON_STATICFORLOOP_HH

#include <cstddef>
#include <type_traits>

namespace Dune
{
  // implementation of staticForLoop
  namespace Impl
  {
    template<class ST, ST begin, ST end, bool condition>
    struct StaticForLoopImpl;

    template<class ST, ST begin, ST end>
    using StaticForLoop
      = StaticForLoopImpl<ST, begin, end, (begin < end)>;

    template<class ST, ST begin, ST end>
    struct StaticForLoopImpl<ST, begin, end, true>
    {
      using Next = StaticForLoop<ST, begin+1, end>;

      template<class F, class ... Args>
      static void apply(F&& f, Args&&... args)
      {
        f(std::integral_constant<ST, begin>(), std::forward<Args>(args) ...);
        Next::apply(std::forward<F>(f), std::forward<Args>(args) ...);
      }
    };

    template<class ST, ST begin, ST end>
    struct StaticForLoopImpl<ST, begin, end, false>
    {
      template<class F, class ... Args>
      static void apply(F&& f, Args&&...)
      {}
    };

  } //end namespace Impl



  /**
   * \brief Static for loop
   *
   * \ingroup Utility
   *
   * Run static for-loop from 'begin' to 'end-1' with functor.
   * The functor is called with \p std::index_constant<i>
   * as first argument. All other arguments of this method
   * are forwarded to the functor.
   */
  template<std::size_t begin_t, std::size_t end_t, class F, class ... Args>
  void staticForLoop(F&& f, Args&&... args)
  {
    using Runner = Impl::StaticForLoop<std::size_t, begin_t, end_t>;
    Runner::apply(std::forward<F>(f), std::forward<Args>(args) ...);
  }


  // implementation of staticFindInRange
  namespace Impl
  {
    template<class ST, ST begin, ST end, bool condition>
    struct StaticFindInRangeImpl;

    template<class ST, ST begin, ST end>
    using StaticFindInRange
      = StaticFindInRangeImpl<ST, begin, end, (begin < end)>;

    template<class ST, ST begin, ST end>
    struct StaticFindInRangeImpl<ST, begin, end, true>
    {
      using Next = StaticFindInRange<ST, begin+1, end>;

      template<class F, class ... Args>
      static ST apply(F&& f, Args&&... args)
      {
        if (f(std::integral_constant<ST, begin>(), std::forward<Args>(args) ...))
          return begin;
        return Next::apply(std::forward<F>(f), std::forward<Args>(args) ...);
      }
    };

    template<class ST, ST begin, ST end>
    struct StaticFindInRangeImpl<ST, begin, end, false>
    {
      template<class F, class ... Args>
      static ST apply(F&& f, Args&&...)
      {
        return end;
      }
    };

  } //end namespace Impl



  /**
   * \brief Static find loop
   *
   * \ingroup Utility
   *
   * Run static for-loop from 'begin' to 'end-1' with functor.
   * The functor is called with \p std::index_constant<i>
   * as first argument. All other arguments of this method
   * are forwarded to the functor. If the functor returns
   * true the loop is terminated. The loop returns the index of
   * termination.
   */
  template<std::size_t begin_t, std::size_t end_t, class F, class ... Args>
  std::size_t staticFindInRange(F&& f, Args&&... args)
  {
    using Runner = Impl::StaticFindInRange<std::size_t, begin_t, end_t>;
    return Runner::apply(std::forward<F>(f), std::forward<Args>(args) ...);
  }


  // implementation if staticAccumulate
  namespace Impl
  {
    template<class ST, ST begin, ST end, bool condition>
    struct StaticAccumulateImpl;

    template<class ST, ST begin, ST end>
    using StaticAccumulate
      = StaticAccumulateImpl<ST, begin, end, (begin < end)>;

    template<class ST, ST begin, ST end>
    struct StaticAccumulateImpl<ST, begin, end, true>
    {
      using Next = StaticAccumulate<ST, begin+1, end>;

      template<class F, class T, class BinaryOp, class ... Args>
      static auto eval(F&& f, T init, BinaryOp&& op, Args&&... args)
      {
        return op( f(std::integral_constant<ST, begin>(), std::forward<Args>(args)...),
                   Next::eval(std::forward<F>(f), init,
                              std::forward<BinaryOp>(op),
                              std::forward<Args>(args)...) );
      }
    };

    template<class ST, ST begin, ST end>
    struct StaticAccumulateImpl<ST, begin, end, false>
    {
      template<class F, class T, class BinaryOp, class ... Args>
      static auto eval(F&& f, T init, BinaryOp&& op, Args&&... args)
      {
        return init;
      }
    };

  } //end namespace Impl



  /**
   * \brief Static accumulate
   *
   * \ingroup Utility
   *
   * Run static accumulate-loop from 'begin' to 'end-1' using an binary operation
   * op applied to the functor values, when called with \p std::index_constant<i>.
   * The break value for the accumulation is the value \param init.
   */
  template<std::size_t begin_t, std::size_t end_t, class F, class T, class BinaryOp, class ... Args>
  auto staticAccumulate(F&& f, T init, BinaryOp&& op, Args&&... args)
  {
    using Runner = Impl::StaticAccumulate<std::size_t, begin_t, end_t>;
    return Runner::eval(std::forward<F>(f), init,
                        std::forward<BinaryOp>(op),
                        std::forward<Args>(args) ...);
  }
} // namespace Dune



#endif //DUNE_LOCALFUNCTIONS_COMMON_STATICFORLOOP_HH
