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
    namespace Impl
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
          lb.size(),
          requireConvertible<unsigned int>( lb.size() ),

          //! \brief Polynomial order of the shape functions
          lb.order(),
          requireConvertible<unsigned int>( lb.order() ),

          //! \brief Evaluate all basis function at given position.
          lb.evaluateFunction(std::declval<D>(), vectorR),

          //! \brief Evaluate jacobian of all shape functions at given position.
          lb.evaluateJacobian(std::declval<D>(), vectorJ),

          //! \brief Evaluate partial derivatives of any order of all shape
          //! functions, using a multiIndex notation.
          lb.partial(std::array<unsigned int, n>{}, std::declval<D>(), vectorR),

          //! \brief Evaluate partial derivatives of any order of all shape
          //! functions, using a directions-array. \deprecated
          lb.template evaluate<dOrder>(std::array<int, dOrder>{}, std::declval<D>(), vectorR),
        0);
      };
    } // end namespace Impl


    struct LocalBasisConcept
    {
      template<class LB>
      auto require(LB&&) -> decltype
      (
        requireType<typename LB::Traits>(),
        requireConcept<Impl::LocalBasisConcept<typename LB::Traits>, LB>()
      );
    };

    //! \brief A static concept check that returns true, if the type LB models a
    //! \ref LocalBasisConcept.
    template <class LB>
    static constexpr bool isLocalBasis()
    {
      return Dune::models< LocalBasisConcept, LB >();
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

    //! \brief A static concept check that returns true, if the type LB has a
    //! method evaluateFunction().
    template <class LB>
    static constexpr bool localBasisHasEvaluateFunction()
    {
      using C = LocalBasisPartialConcept<typename LB::Traits>;
      return Dune::models< typename C::EvaluateFunction, LB >();
    }

    //! \brief A static concept check that returns true, if the type LB has a
    //! method evaluateJacobian().
    template <class LB>
    static constexpr bool localBasisHasEvaluateJacobian()
    {
      using C = LocalBasisPartialConcept<typename LB::Traits>;
      return Dune::models< typename C::EvaluateJacobian, LB >();
    }

    //! \brief A static concept check that returns true, if the type LB has a
    //! method partial().
    template <class LB>
    static constexpr bool localBasisHasPartial()
    {
      using C = LocalBasisPartialConcept<typename LB::Traits>;
      return Dune::models< typename C::Partial, LB >();
    }

    //! \brief A static concept check that returns true, if the type LB has a
    //! method evaluate().
    template <class LB>
    static constexpr bool localBasisHasEvaluate()
    {
      using C = LocalBasisPartialConcept<typename LB::Traits>;
      return Dune::models< typename C::Evaluate, LB >();
    }


    // -------------------------------------------------------------------------


    template<class DomainType, class RangeType>
    struct LocalInterpolationConcept
    {
      //! type of virtual function to interpolate
      using FunctionType = Dune::VirtualFunction<DomainType, RangeType>;

      //! type of the coefficient vector in the interpolate method
      using CoefficientType = typename RangeType::field_type;

      std::vector<CoefficientType> out;

      template<class LI>
      auto require(LI&& li) -> decltype
      (
        //! \brief determine coefficients interpolating a given function
        li.interpolate (std::declval<FunctionType>(), out)
      );

    };

    //! \brief A static concept check that returns true, if the type LI models a
    //! \ref LocalInterpolationConcept.
    template <class LI,
              class DomainType = typename LI::FunctionType::DomainType,
              class RangeType = typename LI::FunctionType::RangeType>
    static constexpr bool isLocalInterpolation()
    {
      return Dune::models< LocalInterpolationConcept<DomainType, RangeType>, LI >();
    }


    // -------------------------------------------------------------------------


    struct LocalCoefficientsConcept
    {
      template<class LC>
      auto require(LC&& lc) -> decltype
      (
        //! \brief number of coefficients
        lc.size(), // -> std::size_t
        requireConvertible<std::size_t>( lc.size() ),

        //! \brief get i'th index
        lc.localKey(std::declval<std::size_t>()),
      0);

    };

    //! \brief A static concept check that returns true, if the type LC models a
    //! \ref LocalCoefficientsConcept.
    template <class LC>
    static constexpr bool isLocalCoefficients()
    {
      return Dune::models< LocalCoefficientsConcept, LC >();
    }


    // -------------------------------------------------------------------------


    namespace Impl
    {
      template<class Traits>
      struct LocalFiniteElementConcept;

      template<class LB, class LC, class LI>
      struct LocalFiniteElementConcept< LocalFiniteElementTraits<LB, LC, LI> >
      {
        using DomainType = typename LB::Traits::DomainType;
        using RangeType = typename LB::Traits::RangeType;

        template<class FE>
        auto require(FE&& fe) -> decltype
        (
          //! \brief Return `Traits::LocalBasisType const&`
          fe.localBasis(),
          requireConcept<LocalBasisConcept<typename LB::Traits>>( fe.localBasis() ),

          //! \brief Return `Traits::LocalCoefficientsType const&`
          fe.localCoefficients(),
          requireConcept<LocalCoefficientsConcept>( fe.localCoefficients() ),

          //! \brief Return `Traits::LocalInterpolationType const&`
          fe.localInterpolation(),
          requireConcept<LocalInterpolationConcept<DomainType, RangeType>>( fe.localInterpolation() ),

          //! \brief Number of shape functions in this finite element
          fe.size(),
          requireConvertible<unsigned int>( fe.size() ),

          //! \brief Return `const GeometryType`
          fe.type(),
          requireConvertible<GeometryType>( fe.type() ),
        0);

      };
    } // end namespace Impl


    struct LocalFiniteElementConcept
    {
      template<class LFE>
      auto require(LFE&&) -> decltype
      (
        requireType<typename LFE::Traits>(),
        requireConcept<Impl::LocalFiniteElementConcept<typename LFE::Traits>, LFE>()
      );
    };

    //! \brief A static concept check that returns true, if the type LFE models a
    //! \ref LocalFiniteElementConcept.
    template <class LFE>
    static constexpr bool isLocalFiniteElement()
    {
      return Dune::models< LocalFiniteElementConcept, LFE >();
    }

  } // end namespace Concept
} // namespace Dune

#endif //DUNE_LOCALFUNCTIONS_COMMON_CONCEPTS_HH
