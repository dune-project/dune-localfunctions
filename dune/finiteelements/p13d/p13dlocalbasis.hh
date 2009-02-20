// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_P13DLOCALBASIS_HH
#define DUNE_P13DLOCALBASIS_HH

#include "../common/localbasis.hh"

namespace Dune
{
  /**@ingroup LocalBasisImplementation
         \brief Linear Lagrange shape functions on the tetrahedron.

         Defines the linear shape functions on tetrahedron.

         \tparam D Type to represent the field in the domain.
         \tparam R Type to represent the field in the range.

         \nosubgrouping
   */
  template<class D, class R>
  class P13DLocalBasis :
    public C1LocalBasisInterface<
        C1LocalBasisTraits<D,3,Dune::FieldVector<D,3>,R,1,Dune::FieldVector<R,1>,
            Dune::FieldVector<Dune::FieldVector<R,3>,1> >,
        P13DLocalBasis<D,R> >
  {
  public:
    //! \brief export type traits for function signature
    typedef C1LocalBasisTraits<D,3,Dune::FieldVector<D,3>,R,1,Dune::FieldVector<R,1>,
        Dune::FieldVector<Dune::FieldVector<R,3>,1> > Traits;

    //! \brief number of shape functions
    unsigned int size () const
    {
      return 4;
    }

    //! \brief Evaluate all shape functions
    inline void evaluateFunction (const typename Traits::DomainType& in,
                                  std::vector<typename Traits::RangeType>& out) const
    {
      out.resize(4);
      out[0] = 1.0-in[0]-in[1]-in[2];
      out[1] = in[0];
      out[2] = in[1];
      out[3] = in[2];
    }

    //! \brief Evaluate Jacobian of all shape functions
    inline void
    evaluateJacobian (const typename Traits::DomainType& in,         // position
                      std::vector<typename Traits::JacobianType>& out) const      // return value
    {
      out.resize(4);
      out[0][0][0] = -1; out[0][0][1] = -1; out[0][0][2] = -1; // basis function 0
      out[1][0][0] =  1; out[1][0][1] =  0; out[1][0][2] =  0; // basis function 1
      out[2][0][0] =  0; out[2][0][1] =  1; out[2][0][2] =  0; // basis function 2
      out[3][0][0] =  0; out[3][0][1] =  0; out[3][0][2] =  1; // basis function 3
    }

    //! \brief Polynomial order of the shape functions
    unsigned int order () const
    {
      return 1;
    }
  };
}
#endif
