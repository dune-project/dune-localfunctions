// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LOCALFUNCTIONS_COMMON_STATICFORLOOP_HH
#define DUNE_LOCALFUNCTIONS_COMMON_STATICFORLOOP_HH

#include <cstddef>
#include <type_trits>

namespace Dune
{
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
      template<class F, class ... Args>
      static void apply(F&& f, Args&&... args)
      {
        f(std::integral_constant<ST, begin>(), std::forward<Args>(args) ...);
        StaticForLoop<ST, begin+1, end>::apply(std::forward<F>(f), std::forward<Args>(args) ...);
      }
    };

    template<class ST, ST begin, ST end>
    struct StaticForLoop<ST, begin, end, false>
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
    Impl::StaticForLoop<std::size_t, begin_t, end_t>::apply(std::forward<F>(f), std::forward<Args>(args) ...);
  }



  namespace Impl
  {
    template<class ST, ST begin, ST end, bool condition>
    struct StaticFindInRange;

    template<class ST, ST begin, ST end>
    using StaticFindInRange
            = StaticFindInRangeImpl<ST, begin, end, (begin < end)>;

    template<class ST, ST begin, ST end>
    struct StaticFindInRangeImpl<ST, begin, end, true>
    {
      template<class F, class ... Args>
      static ST apply(F&& f, Args&&... args)
      {
        if (f(std::integral_constant<ST, begin>(), std::forward<Args>(args) ...))
          return begin;
        return StaticFindInRange<ST, begin+1, end>::apply(std::forward<F>(f), std::forward<Args>(args) ...);
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
   * true the loop is terminated. The loop returns true when
   * easy terminating and false otherwise.
   */
  template<std::size_t begin_t, std::size_t end_t, class F, class ... Args>
  std::size_t staticFindInRange(F&& f, Args&&... args)
  {
    return Impl::StaticFindInRange<std::size_t, begin_t, end_t>::apply(std::forward<F>(f), std::forward<Args>(args) ...);
  }


} // namespace Dune



#endif //DUNE_LOCALFUNCTIONS_COMMON_STATICFORLOOP_HH
