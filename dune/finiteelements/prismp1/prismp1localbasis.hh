// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PRISM_P1_LOCALBASIS_HH
#define DUNE_PRISM_P1_LOCALBASIS_HH

#include "../common/localbasis.hh"

namespace Dune
{
  /**@ingroup LocalBasisImplementation
         \brief Linear Lagrange shape functions on the prism.

         Defines the linear shape functions on the prism.

         \tparam D Type to represent the field in the domain.
         \tparam R Type to represent the field in the range.

         \nosubgrouping
   */
  template<class D, class R>
  class PrismP1LocalBasis :
    public C1LocalBasisInterface<
        C1LocalBasisTraits<D,3,Dune::FieldVector<D,3>,R,1,Dune::FieldVector<R,1>,
            Dune::FieldVector<Dune::FieldVector<R,3>,1> >
#ifndef DUNE_VIRTUAL_SHAPEFUNCTIONS
        , PrismP1LocalBasis<D,R>
#endif
        >
  {
  public:
    //! \brief export type traits for function signature
    typedef C1LocalBasisTraits<D,3,Dune::FieldVector<D,3>,R,1,Dune::FieldVector<R,1>,
        Dune::FieldVector<Dune::FieldVector<R,3>,1> > Traits;

    //! \brief number of shape functions
    unsigned int size () const
    {
      return 6;
    }

    //! \brief Evaluate all shape functions
    inline void evaluateFunction (const typename Traits::DomainType& in,
                                  std::vector<typename Traits::RangeType>& out) const
    {
      out.resize(6);
      out[0] = (1.0-in[0]-in[1])*(1.0-in[2]);
      out[1] = in[0]*(1-in[2]);
      out[2] = in[1]*(1-in[2]);
      out[3] = in[2]*(1.0-in[0]-in[1]);
      out[4] = in[0]*in[2];
      out[5] = in[1]*in[2];
    }

    //! \brief Evaluate Jacobian of all shape functions
    inline void
    evaluateJacobian (const typename Traits::DomainType& in,         // position
                      std::vector<typename Traits::JacobianType>& out) const      // return value
    {
      out.resize(6);
      out[0][0][0] = in[2]-1; out[0][0][1] = in[2]-1; out[0][0][2] = in[0]+in[1]-1;  // basis function 0
      out[1][0][0] = 1-in[2]; out[1][0][1] = 0;       out[1][0][2] = -in[0];         // basis function 1
      out[2][0][0] = 0;       out[2][0][1] = 1-in[2]; out[2][0][2] = -in[1];         // basis function 2
      out[3][0][0] = -in[2];  out[3][0][1] = -in[2];  out[3][0][2] = 1-in[0]-in[1];  // basis function 3
      out[4][0][0] = in[2];   out[4][0][1] = 0;       out[4][0][2] = in[0];          // basis function 4
      out[5][0][0] = 0;       out[5][0][1] = in[2];   out[5][0][2] = in[1];          // basis function 5
    }

    //! \brief Polynomial order of the shape functions
    unsigned int order () const
    {
      return 1;
    }
  };
}
#endif
