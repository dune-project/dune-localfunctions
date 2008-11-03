// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_P0LOCALBASIS_HH
#define DUNE_P0LOCALBASIS_HH

#include "../common/localbasis.hh"

namespace Dune
{
  /**@ingroup LocalBasisImplementation
         \brief Constant shape function

         Defines the constant scalar shape function in d dimensions. Is
         valid on any type of reference element.

         - <tt>D</tt>: Type to represent the field in the domain.
         - <tt>R</tt>: Type to represent the field in the range.
         - <tt>d</tt>: Domain dimension

         \nosubgrouping
   */
  template<class D, class R, int d>
  class P0LocalBasis :
    public C1LocalBasisInterface<
        C1LocalBasisTraits<D,d,Dune::FieldVector<D,d>,R,1,Dune::FieldVector<R,1>,
            Dune::FieldVector<Dune::FieldVector<R,d>,1> >,
        P0LocalBasis<D,R,d> >
  {
  public:
    //! \brief export type traits for function signature
    typedef C1LocalBasisTraits<D,d,Dune::FieldVector<D,d>,R,1,Dune::FieldVector<R,1>,
        Dune::FieldVector<Dune::FieldVector<R,d>,1> > Traits;

    //! \brief number of shape functions
    unsigned int size () const
    {
      return 1;
    }

    //! \brief Evaluate all shape functions
    inline void evaluateFunction (const typename Traits::DomainType& in,
                                  std::vector<typename Traits::RangeType>& out) const
    {
      out.resize(1);
      out[0] = 1;
    }

    //! \brief Evaluate Jacobian of all shape functions
    inline void
    evaluateJacobian (const typename Traits::DomainType& in,       // position
                      std::vector<typename Traits::JacobianType>& out) const                        // return value
    {
      out.resize(1);
      for (int i=0; i<d; i++)
        out[0][0][i] = 0;
    }

    //! \brief Local interpolation of a function
    template<typename E, typename F, typename C>
    void interpolate (const E& e, const F& f, std::vector<C>& out) const
    {
      typedef typename Traits::DomainType DomainType;
      typedef typename Traits::RangeType RangeType;
      typedef typename Traits::DomainFieldType DF;
      const int dim=Traits::dimDomain;

      DomainType x = Dune::ReferenceElements<DF,dim>::general(e.type()).position(0,0);
      RangeType y;

      out.resize(1);
      f.eval_local(e,x,y); out[0] = y;
    }

    //! \brief Polynomial order of the shape functions
    unsigned int order () const
    {
      return 0;
    }
  };

}

#endif
