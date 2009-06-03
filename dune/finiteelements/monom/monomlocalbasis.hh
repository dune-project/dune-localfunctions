// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_MONOMLOCALBASIS_HH
#define DUNE_MONOMLOCALBASIS_HH

#include <dune/grid/common/referenceelements.hh>

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

    //! Access output vector of evaluateFunction() and evaluate()
    template <typename Traits>
    class EvalAccess {
      std::vector<typename Traits::RangeType> &out;
    public:
      EvalAccess(std::vector<typename Traits::RangeType> &out_)
        : out(out_)
      { }
      typename Traits::RangeFieldType &operator[](unsigned int index)
      {
        return out[index][0];
      }
    };

    //! Access output vector of evaluateJacobian()
    template <typename Traits>
    class JacobianAccess {
      std::vector<typename Traits::JacobianType> &out;
      unsigned int row;
    public:
      JacobianAccess(std::vector<typename Traits::JacobianType> &out_,
                     unsigned int row_)
        : out(out_), row(row_)
      { }
      typename Traits::RangeFieldType &operator[](unsigned int index)
      {
        return out[index][row][0];
      }
    };

    /** Template Metaprogramm for evaluating monomial shapefunctions
     *  \internal
     */
    template <typename Traits, int c>
    struct Evaluate
    {
      enum { d = Traits::dimDomain - c };
      template <typename Access>
      static void eval (const typename Traits::DomainType &in,
                        const array<int, Traits::dimDomain> &derivatives,
                        typename Traits::RangeFieldType prod,
                        int bound, int& index, Access &access)
      {
        for (int newbound=0; newbound<=bound; newbound++)
        {
          int e = bound-newbound;
          /*if(e < derivatives[d])
             Evaluate<Traits,c-1>::
              eval(in, derivatives, 0, newbound, index, access);
             else*/{
            /*int divisor = 1;
               for(int i = e - derivatives[d] + 1; i <= e; ++i)
               divisor *= i;*/
            Evaluate<Traits,c-1>::
            eval(in, derivatives,
                 prod*std::pow(in[d], typename Traits::DomainFieldType(e /*-derivatives[d]*/)) /*/ divisor* /,
                 newbound, index, access);
          }
        }
      }
    };

    /** \copydoc evaluate
     *  \brief Specializes the end of the recursion
     *  \internal
     */
    template <typename Traits>
    struct Evaluate<Traits, 1>
    {
      enum { d = Traits::dimDomain-1 };
      template <typename Access>
      static void eval (const typename Traits::DomainType &in,
                        const array<int, Traits::dimDomain> &derivatives,
                        typename Traits::RangeFieldType prod,
                        int bound, int& index, Access &access)
      {
        /*if(bound < derivatives[d])
           prod = 0;
           else*/{
          /*int divisor = 1;
             for(int i = bound - derivatives[d] + 1; i <= bound; ++i)
             divisor *= i;*/
          prod *= std::pow(in[d], typename Traits::DomainFieldType(bound /*-derivatives[d]*/)) /*/ divisor* /;
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

         - <tt>D</tt>: Type to represent the field in the domain.
         - <tt>R</tt>: Type to represent the field in the range.
         - <tt>d</tt>: Domain dimension
     - <tt>p</tt>: polynomial order of the shapefunctions

         \nosubgrouping
   */
  template<class D, class R, unsigned int d, unsigned int p>
  class MonomLocalBasis :
    public CkLocalBasisInterface<
        CkLocalBasisTraits<D,d,Dune::FieldVector<D,d>,
            R,1,Dune::FieldVector<R,1>,
            Dune::FieldVector<Dune::FieldVector<R,d>,1>, p>,
        MonomLocalBasis<D,R,d,p> >
  {
    enum { static_size = MonomImp::Size<d,p>::val };

  public:
    //! \brief export type traits for function signature
    typedef CkLocalBasisTraits<D,d,Dune::FieldVector<D,d>,R,1,Dune::FieldVector<R,1>,
        Dune::FieldVector<Dune::FieldVector<R,d>,1>,p> Traits;

    //! \brief number of shape functions
    unsigned int size () const
    {
      return static_size;
    }

    //! \brief Evaluate all shape functions
    inline void evaluateFunction (const typename Traits::DomainType& in,
                                  std::vector<typename Traits::RangeType>& out) const
    {
      evaluate<0>(array<int, 0>(), in, out);
    }

    //! return given derivative of all components
    template<int k>
    inline void evaluate (const array<int,k>& directions,
                          const typename Traits::DomainType& in,
                          std::vector<typename Traits::RangeType>& out) const
    {
      out.resize(size());
      int index = 0;
      array<int, d> derivatives;
      for(unsigned int i = 0; i < d; ++i) derivatives[i] = 0;
      for(int i = 0; i < k; ++i) ++derivatives[directions[i]];
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
      array<int, d> derivatives;
      for(int i = 0; i < d; ++i)
        derivatives[i] = 0;
      for(int i = 0; i < d; ++i)
      {
        derivatives[i] = 1;
        int index = 0;
        MonomImp::JacobianAccess<Traits> access(out, i);
        MonomImp::Evaluate<Traits, d>::eval(in, derivatives, 1, p, index, access);
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

#endif // DUNE_MONOMLOCALBASIS_HH
