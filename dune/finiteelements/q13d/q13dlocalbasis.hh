// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_Q13DLOCALBASIS_HH
#define DUNE_Q13DLOCALBASIS_HH

#include "../common/localbasis.hh"

namespace Dune
{
  /**@ingroup LocalBasisImplementation
         \brief Lagrange shape functions of order 1 on the reference hexahedron.

         Also known as \f$Q^1\f$.

         - <tt>D</tt>: Type to represent the field in the domain.
         - <tt>R</tt>: Type to represent the field in the range.

         \nosubgrouping
   */
  template<class D, class R>
  class Q13DLocalBasis :
    public C1LocalBasisInterface<
        C1LocalBasisTraits<D,3,Dune::FieldVector<D,3>,R,1,Dune::FieldVector<R,1>,
            Dune::FieldVector<Dune::FieldVector<R,3>,1> >,
        Q13DLocalBasis<D,R> >
  {
  public:
    typedef C1LocalBasisTraits<D,3,Dune::FieldVector<D,3>,R,1,Dune::FieldVector<R,1>,
        Dune::FieldVector<Dune::FieldVector<R,3>,1> > Traits;

    //! \brief number of shape functions
    unsigned int size () const
    {
      return 8;
    }

    //! \brief Evaluate all shape functions
    inline void evaluateFunction (const typename Traits::DomainType& in,
                                  std::vector<typename Traits::RangeType>& out) const
    {
      out.resize(8);
      out[0] = (1-in[0])*(1-in[1])*(1-in[2]);
      out[1] = (  in[0])*(1-in[1])*(1-in[2]);
      out[2] = (1-in[0])*(  in[1])*(1-in[2]);
      out[3] = (  in[0])*(  in[1])*(1-in[2]);
      out[4] = (1-in[0])*(1-in[1])*(  in[2]);
      out[5] = (  in[0])*(1-in[1])*(  in[2]);
      out[6] = (1-in[0])*(  in[1])*(  in[2]);
      out[7] = (  in[0])*(  in[1])*(  in[2]);
    }

    //! \brief Evaluate Jacobian of all shape functions
    inline void
    evaluateJacobian (const typename Traits::DomainType& in,       // position
                      std::vector<typename Traits::JacobianType>& out) const                        // return value
    {
      out.resize(8);
      out[0][0][0] = -(1-in[1])*(1-in[2]); out[0][0][1] = -(1-in[0])*(1-in[2]); out[0][0][2] = -(1-in[0])*(1-in[1]);
      out[1][0][0] = +(1-in[1])*(1-in[2]); out[1][0][1] = -(  in[0])*(1-in[2]); out[1][0][2] = -(  in[0])*(1-in[1]);
      out[2][0][0] = -(  in[1])*(1-in[2]); out[2][0][1] = +(1-in[0])*(1-in[2]); out[2][0][2] = -(1-in[0])*(  in[1]);
      out[3][0][0] = +(  in[1])*(1-in[2]); out[3][0][1] = +(  in[0])*(1-in[2]); out[3][0][2] = -(  in[0])*(  in[1]);
      out[4][0][0] = -(1-in[1])*(  in[2]); out[4][0][1] = -(1-in[0])*(  in[2]); out[4][0][2] = +(1-in[0])*(1-in[1]);
      out[5][0][0] = +(1-in[1])*(  in[2]); out[5][0][1] = -(  in[0])*(  in[2]); out[5][0][2] = +(  in[0])*(1-in[1]);
      out[6][0][0] = -(  in[1])*(  in[2]); out[6][0][1] = +(1-in[0])*(  in[2]); out[6][0][2] = +(1-in[0])*(  in[1]);
      out[7][0][0] = +(  in[1])*(  in[2]); out[7][0][1] = +(  in[0])*(  in[2]); out[7][0][2] = +(  in[0])*(  in[1]);

      //        out[0][0][0] = -(1-in[1])*(1-in[2]); out[0][0][1] = -(1-in[0])*(1-in[2]); out[0][0][2] = -(1-in[1])*(1-in[2]);
      //        out[1][0][0] = +(1-in[1])*(1-in[2]); out[1][0][1] = -(1-in[0])*(1-in[2]); out[1][0][2] = -(1-in[1])*(1-in[2]);
      //        out[2][0][0] = -(1-in[1])*(1-in[2]); out[2][0][1] = +(1-in[0])*(1-in[2]); out[2][0][2] = -(1-in[1])*(1-in[2]);
      //        out[3][0][0] = +(1-in[1])*(1-in[2]); out[3][0][1] = +(1-in[0])*(1-in[2]); out[3][0][2] = -(1-in[1])*(1-in[2]);
      //        out[4][0][0] = -(1-in[1])*(1-in[2]); out[4][0][1] = -(1-in[0])*(1-in[2]); out[4][0][2] = +(1-in[1])*(1-in[2]);
      //        out[5][0][0] = +(1-in[1])*(1-in[2]); out[5][0][1] = -(1-in[0])*(1-in[2]); out[5][0][2] = +(1-in[1])*(1-in[2]);
      //        out[6][0][0] = -(1-in[1])*(1-in[2]); out[6][0][1] = +(1-in[0])*(1-in[2]); out[6][0][2] = +(1-in[1])*(1-in[2]);
      //        out[7][0][0] = +(1-in[1])*(1-in[2]); out[7][0][1] = +(1-in[0])*(1-in[2]); out[7][0][2] = +(1-in[1])*(1-in[2]);
    }

    //! \brief Polynomial order of the shape functions
    unsigned int order () const
    {
      return 1;
    }
  };
}
#endif
