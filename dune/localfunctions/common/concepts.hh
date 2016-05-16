// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LOCALFUNCTIONS_COMMON_CONCEPTS_HH
#define DUNE_LOCALFUNCTIONS_COMMON_CONCEPTS_HH

#include <array>
#include <vector>

#include <dune/common/concept.hh>

namespace Dune
{
  namespace Concept
  {
    template<class Traits>
    struct LocalBasisConcept;

    template<class DF, int n, class D, class RF, int m, class R, class J, int dOrder>
    struct LocalBasisConcept< LocalBasisTraits<DF,n,D,RF,m,R,J,dOrder> >
    {
      std::vector<R> vectorR; // RangeType
      std::vector<J> vectorJ; // JacobianType

      template<class LB>
      auto require(LB&& lb) -> decltype
      (
        //! \brief Number of shape functions
        lb.size(),  // -> unsigned int

        //! \brief Polynomial order of the shape functions
        lb.order(), // -> unsigned int

        //! \brief Evaluate all basis function at given position.
        lb.evaluateFunction(std::declval<D>(), vectorR),

        //! \brief Evaluate jacobian of all shape functions at given position.
        lb.evaluateJacobian(std::declval<D>(), vectorJ),

        //! \brief Evaluate partial derivatives of any order of all shape
        //! functions, using a multiIndex notation.
        lb.partial(std::array<unsigned int, n>{}, std::declval<D>(), vectorR),

        //! \brief Evaluate partial derivatives of any order of all shape
        //! functions, using a directions-array.
        lb.template evaluate<dOrder>(std::array<int, dOrder>{}, std::declval<D>(), vectorR)
      );

    };

    //! \brief A static concept check that returns true, if the type \param LB models a
    //! \ref LocalBasisConcept.
    template <class LB>
    static constexpr bool isLocalBasis()
    {
      return Dune::models< LocalBasisConcept<typename LB::Traits>, LB >();
    }


    template<class Traits>
    struct LocalBasisPartialConcept;

    template<class DF, int n, class D, class RF, int m, class R, class J, int dOrder>
    struct LocalBasisPartialConcept< LocalBasisTraits<DF,n,D,RF,m,R,J,dOrder> >
    {
      // check whether a method evaluateFunction() exists
      struct EvaluateFunction
      {
        std::vector<R> out; // RangeType
        template<class LB>
        auto require(LB&& lb) -> decltype(
          lb.evaluateFunction(std::declval<D>(), out)
        );
      };

      // check whether a method evaluateJacobian() exists
      struct EvaluateJacobian
      {
        std::vector<J> out; // JacobianType
        template<class LB>
        auto require(LB&& lb) -> decltype(
          lb.evaluateJacobian(std::declval<D>(), out)
        );
      };

      // check whether a method partial() exists
      struct Partial
      {
        std::vector<R> out; // RangeType
        template<class LB>
        auto require(LB&& lb) -> decltype(
          lb.partial(std::array<unsigned int, n>{}, std::declval<D>(), out)
        );
      };

      // check whether a method evaluate() exists
      struct Evaluate
      {
        std::vector<R> out; // RangeType
        template<class LB>
        auto require(LB&& lb) -> decltype(
          lb.template evaluate<dOrder>(std::array<int, dOrder>{}, std::declval<D>(), out)
        );
      };
    };

    //! \brief A static concept check that returns true, if the type \param LB has a
    //! method evaluateFunction().
    template <class LB>
    static constexpr bool localBasisHasEvaluateFunction()
    {
      using C = LocalBasisPartialConcept<typename LB::Traits>;
      return Dune::models< typename C::EvaluateFunction, LB >();
    }

    //! \brief A static concept check that returns true, if the type \param LB has a
    //! method evaluateJacobian().
    template <class LB>
    static constexpr bool localBasisHasEvaluateJacobian()
    {
      using C = LocalBasisPartialConcept<typename LB::Traits>;
      return Dune::models< typename C::EvaluateJacobian, LB >();
    }

    //! \brief A static concept check that returns true, if the type \param LB has a
    //! method partial().
    template <class LB>
    static constexpr bool localBasisHasPartial()
    {
      using C = LocalBasisPartialConcept<typename LB::Traits>;
      return Dune::models< typename C::Partial, LB >();
    }

    //! \brief A static concept check that returns true, if the type \param LB has a
    //! method evaluate().
    template <class LB>
    static constexpr bool localBasisHasEvaluate()
    {
      using C = LocalBasisPartialConcept<typename LB::Traits>;
      return Dune::models< typename C::Evaluate, LB >();
    }

  } // end namespace Concept
} // namespace Dune

#endif //DUNE_LOCALFUNCTIONS_COMMON_CONCEPTS_HH
