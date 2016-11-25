// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LOCALFUNCTIONS_MONOMIAL_MONOMIALLOCALBASIS_HH
#define DUNE_LOCALFUNCTIONS_MONOMIAL_MONOMIALLOCALBASIS_HH

#include <array>
#include <cassert>
#include <numeric>

#include <dune/common/fmatrix.hh>
#include <dune/common/deprecated.hh>

#include "../common/localbasis.hh"

namespace Dune
{
  namespace MonomImp {
    /** template meta program to calculate the number of shape functions
     *  \internal
     */
    template<int d, int k>
    struct Size {
      enum { val = Size<d,k-1>::val+Size<d-1,k>::val };
    };
    template<int d>
    struct Size<d, 0> {
      enum { val = 1 };
    };
    template<int k>
    struct Size<0, k> {
      enum { val = 1 };
    };
    template<>
    struct Size<0, 0> {
      enum { val = 1 };
    };

    template<class T>
    T ipow(T base, int exp)
    {
      T result(1);
      while (exp)
      {
        if (exp & 1)
          result *= base;
        exp >>= 1;
        base *= base;
      }
      return result;
    }

    //! Access output vector of evaluateFunction() and evaluate()
    template <typename Traits>
    class EvalAccess {
      std::vector<typename Traits::RangeType> &out;
#ifndef NDEBUG
      unsigned int first_unused_index;
#endif

    public:
      EvalAccess(std::vector<typename Traits::RangeType> &out_)
        : out(out_)
#ifndef NDEBUG
          , first_unused_index(0)
#endif
      { }
#ifndef NDEBUG
      ~EvalAccess() {
        assert(first_unused_index == out.size());
      }
#endif
      typename Traits::RangeFieldType &operator[](unsigned int index)
      {
        assert(index < out.size());
#ifndef NDEBUG
        if(first_unused_index <= index)
          first_unused_index = index+1;
#endif
        return out[index][0];
      }
    };

    //! Access output vector of evaluateJacobian()
    template <typename Traits>
    class JacobianAccess {
      std::vector<typename Traits::JacobianType> &out;
      unsigned int row;
#ifndef NDEBUG
      unsigned int first_unused_index;
#endif

    public:
      JacobianAccess(std::vector<typename Traits::JacobianType> &out_,
                     unsigned int row_)
        : out(out_), row(row_)
#ifndef NDEBUG
          , first_unused_index(0)
#endif
      { }
#ifndef NDEBUG
      ~JacobianAccess() {
        assert(first_unused_index == out.size());
      }
#endif
      typename Traits::RangeFieldType &operator[](unsigned int index)
      {
        assert(index < out.size());
#ifndef NDEBUG
        if(first_unused_index <= index)
          first_unused_index = index+1;
#endif
        return out[index][0][row];
      }
    };

    /** Template Metaprogramm for evaluating monomial shapefunctions
     *  \internal
     *
     *  \tparam Traits The Traits class of the monomial shape functions to
     *                 evaluate -- used to get DomainType etc.
     *  \tparam c      The "codim of the next dimension to try for factors".
     *                 Unfortunately, we cannot recurs over that dimension
     *                 directly, since the end of the recursion cannot be
     *                 specialized for dimDomain-1, but we can recurs over
     *                 dimDomain minus that dimension, since it can be
     *                 specialized for 1.
     */
    template <typename Traits, int c>
    struct Evaluate
    {
      enum {
        //! The next dimension to try for factors
        d = Traits::dimDomain - c
      };
      /** \todo
       *
       *  \tparam Access Wrapper around the result vector, so we don't have to
       *                 copy the output and can still use the same code for
       *                 both the usual drivatives and for the Jacobian
       */
      template <typename Access>
      static void eval ( //! The point at which to evaluate
        const typename Traits::DomainType &in,
        //! The number of partial derivatives, one entry for
        //! each dimension
        const std::array<int, Traits::dimDomain> &derivatives,
        //! The product accumulated for the dimensions which
        //! have already been handled
        typename Traits::RangeFieldType prod,
        //! The number of factors still to go
        int bound,
        //! The index of the next entry in the output to fill
        int& index,
        //! The wrapper used to access the output vector
        Access &access)
      {
        // start with the highest exponent for this dimension, then work down
        for (int e = bound; e >= 0; --e)
        {
          // the rest of the available exponents, to be used by the other
          // dimensions
          int newbound = bound - e;
          if(e < derivatives[d])
            Evaluate<Traits,c-1>::
            eval(in, derivatives, 0, newbound, index, access);
          else {
            int coeff = 1;
            for(int i = e - derivatives[d] + 1; i <= e; ++i)
              coeff *= i;
            // call the evaluator for the next dimension
            Evaluate<Traits,c-1>::
            eval(  // pass the coordinate and the derivatives unchanged
              in, derivatives,
              // also pass the product accumulated so far, but also
              // include the current dimension
              prod * ipow(in[d], e-derivatives[d]) * coeff,
              // pass the number of remaining exponents to the next
              // dimension
              newbound,
              // pass the next index to fill and the output access
              // wrapper
              index, access);
          }
        }
      }
    };

    /** \copydoc Evaluate
     *  \brief Specializes the end of the recursion
     *  \internal
     */
    template <typename Traits>
    struct Evaluate<Traits, 1>
    {
      enum { d = Traits::dimDomain-1 };
      //! \copydoc Evaluate::eval
      template <typename Access>
      static void eval (const typename Traits::DomainType &in,
                        const std::array<int, Traits::dimDomain> &derivatives,
                        typename Traits::RangeFieldType prod,
                        int bound, int& index, Access &access)
      {
        if(bound < derivatives[d])
          prod = 0;
        else {
          int coeff = 1;
          for(int i = bound - derivatives[d] + 1; i <= bound; ++i)
            coeff *= i;
          prod *= ipow(in[d], bound-derivatives[d]) * coeff;
        }
        access[index] = prod;
        ++index;
      }
    };

  } //namespace MonomImp

  /**@ingroup LocalBasisImplementation
         \brief Constant shape function

         Defines the constant scalar shape function in d dimensions. Is
         valid on any type of reference element.

         \tparam D Type to represent the field in the domain.
         \tparam R Type to represent the field in the range.
         \tparam d Domain dimension
     \tparam p polynomial order of the shapefunctions
     \tparam diffOrder Maximum differentiation order to report in the traits.

         \nosubgrouping
   */
  template<class D, class R, unsigned int d, unsigned int p, unsigned diffOrder = p>
  class MonomialLocalBasis
  {
  public:
    //! \brief export type traits for function signature
    typedef LocalBasisTraits<D,d,Dune::FieldVector<D,d>,R,1,Dune::FieldVector<R,1>,
        Dune::FieldMatrix<R,1,d>,diffOrder> Traits;

    //! \brief number of shape functions
    unsigned int size () const
    {
      return MonomImp::Size<d,p>::val;
    }

    //! \brief Evaluate all shape functions
    inline void evaluateFunction (const typename Traits::DomainType& in,
                                  std::vector<typename Traits::RangeType>& out) const
    {
      DUNE_NO_DEPRECATED_BEGIN
      evaluate<0>(std::array<int, 0>(), in, out);
      DUNE_NO_DEPRECATED_END
    }

    /** \brief Evaluate partial derivatives of any order of all shape functions
     * \param order Order of the partial derivatives, in the classic multi-index notation
     * \param in Position where to evaluate the derivatives
     * \param[out] out Return value: the desired partial derivatives
     */
    inline void partial(const std::array<unsigned int,d>& order,
                        const typename Traits::DomainType& in,
                        std::vector<typename Traits::RangeType>& out) const
    {
      auto totalOrder = std::accumulate(order.begin(), order.end(), 0);

      switch (totalOrder)
      {
        case 0:
          evaluateFunction(in,out);
          break;
        case 1:
        {
          std::array<int,1> directions;
          directions[0] = std::find(order.begin(), order.end(), 1)-order.begin();
          DUNE_NO_DEPRECATED_BEGIN
          evaluate<1>(directions, in, out);
          DUNE_NO_DEPRECATED_END
          break;
        }
        case 2:
        {
          std::array<int,2> directions;
          unsigned int counter = 0;
          auto nonconstOrder = order;  // need a copy that I can modify
          for (unsigned int i=0; i<d; i++)
          {
            while (nonconstOrder[i])
            {
              directions[counter++] = i;
              nonconstOrder[i]--;
            }
          }

          DUNE_NO_DEPRECATED_BEGIN
          evaluate<2>(directions, in, out);
          DUNE_NO_DEPRECATED_END
          break;
        }
        default:
          // \todo The 'evaluate' method implements higher derivatives
          DUNE_THROW(NotImplemented, "Desired derivative order is not implemented");
      }
    }

    //! return given derivative of all components
    template<unsigned int k>
    inline void DUNE_DEPRECATED_MSG("Use method 'partial' instead!")
    evaluate (const std::array<int,k>& directions,
                          const typename Traits::DomainType& in,
                          std::vector<typename Traits::RangeType>& out) const
    {
      out.resize(size());
      int index = 0;
      std::array<int, d> derivatives;
      for(unsigned int i = 0; i < d; ++i) derivatives[i] = 0;
      for(unsigned int i = 0; i < k; ++i) ++derivatives[directions[i]];
      MonomImp::EvalAccess<Traits> access(out);
      for(unsigned int lp = 0; lp <= p; ++lp)
        MonomImp::Evaluate<Traits, d>::eval(in, derivatives, 1, lp, index,
                                            access);
    }

    //! \brief Evaluate Jacobian of all shape functions
    inline void
    evaluateJacobian (const typename Traits::DomainType& in,         // position
                      std::vector<typename Traits::JacobianType>& out) const      // return value
    {
      out.resize(size());
      std::array<int, d> derivatives;
      for(unsigned int i = 0; i < d; ++i)
        derivatives[i] = 0;
      for(unsigned int i = 0; i < d; ++i)
      {
        derivatives[i] = 1;
        int index = 0;
        MonomImp::JacobianAccess<Traits> access(out, i);
        for(unsigned int lp = 0; lp <= p; ++lp)
          MonomImp::Evaluate<Traits, d>::eval(in, derivatives, 1, lp, index, access);
        derivatives[i] = 0;
      }
    }

    //! \brief Polynomial order of the shape functions
    unsigned int order () const
    {
      return p;
    }
  };

}

#endif // DUNE_LOCALFUNCTIONS_MONOMIAL_MONOMIALLOCALBASIS_HH
