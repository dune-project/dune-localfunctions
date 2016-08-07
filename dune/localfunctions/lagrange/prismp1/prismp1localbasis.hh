// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PRISM_P1_LOCALBASIS_HH
#define DUNE_PRISM_P1_LOCALBASIS_HH

#include <numeric>

#include <dune/common/fmatrix.hh>

#include <dune/localfunctions/common/localbasis.hh>

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
  class PrismP1LocalBasis
  {
  public:
    //! \brief export type traits for function signature
    typedef LocalBasisTraits<D,3,Dune::FieldVector<D,3>,R,1,Dune::FieldVector<R,1>,
        Dune::FieldMatrix<R,1,3> > Traits;

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

    //! \brief Evaluate partial derivatives of all shape functions
    void partial (const std::array<unsigned int, 3>& order,
                  const typename Traits::DomainType& in,         // position
                  std::vector<typename Traits::RangeType>& out) const      // return value
    {
      auto totalOrder = std::accumulate(order.begin(), order.end(), 0);
      if (totalOrder == 0) {
        evaluateFunction(in, out);
      } else if (totalOrder == 1) {
        auto direction = std::distance(order.begin(), std::find(order.begin(), order.end(), 1));
        out.resize(size());

        switch (direction) {
          case 0:
            out[0] = in[2]-1;
            out[1] = 1-in[2];
            out[2] = 0;
            out[3] = -in[2];
            out[4] = in[2];
            out[5] = 0;
            break;
          case 1:
            out[0] = in[2]-1;
            out[1] = 0;
            out[2] = 1-in[2];
            out[3] = -in[2];
            out[4] = 0;
            out[5] = in[2];
            break;
          case 2:
            out[0] = in[0]+in[1]-1;
            out[1] = -in[0];
            out[2] = -in[1];
            out[3] = 1-in[0]-in[1];
            out[4] = in[0];
            out[5] = in[1];
            break;
          default:
            DUNE_THROW(RangeError, "Component out of range.");
        }
      } else if (totalOrder == 2) {
        out.resize(size());
        if (order[0] == 1 && order[2] == 1) {
          out[0] = 1;
          out[1] =-1;
          out[2] = 0;
          out[3] =-1;
          out[4] = 1;
          out[5] = 0;
        } else if (order[1] == 1 && order[2] == 1) {
          out[0] = 1;
          out[1] = 0;
          out[2] =-1;
          out[3] =-1;
          out[4] = 0;
          out[5] = 1;
        } else {
          for (std::size_t i = 0; i < size(); ++i)
            out[i] = 0;
        }
      } else {
        out.resize(size());
        for (std::size_t i = 0; i < size(); ++i)
          out[i] = 0;
      }
    }

    //! \brief Polynomial order of the shape functions
    unsigned int order () const
    {
      return 1;
    }
  };
}
#endif
