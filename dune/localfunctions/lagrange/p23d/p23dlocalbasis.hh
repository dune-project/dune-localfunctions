// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_P2_3DLOCALBASIS_HH
#define DUNE_P2_3DLOCALBASIS_HH

#include <numeric>

#include <dune/common/fmatrix.hh>

#include <dune/localfunctions/common/localbasis.hh>

namespace Dune
{
  /**@ingroup LocalBasisImplementation
         \brief Quadratic Lagrange shape functions on the tetrahedron.

         Defines the quadratic shape functions on tetrahedron.

         \tparam D Type to represent the field in the domain.
         \tparam R Type to represent the field in the range.

         \nosubgrouping
   */
  template<class D, class R>
  class P23DLocalBasis
  {
  public:
    //! \brief export type traits for function signature
    typedef LocalBasisTraits<D,3,Dune::FieldVector<D,3>,R,1,Dune::FieldVector<R,1>,
        Dune::FieldMatrix<R,1,3> > Traits;

    //! \brief number of shape functions
    unsigned int size () const
    {
      return 10;
    }

    //! \brief Evaluate all shape functions
    inline void evaluateFunction (const typename Traits::DomainType& in,
                                  std::vector<typename Traits::RangeType>& out) const
    {
      out.resize(10);

      int coeff;
      R a[2], b[3], c[3];

      // case 0:
      coeff=2;
      a[0]=1.0;
      a[1]=0.5;
      b[0]=-1.0;
      b[1]=-1.0;
      b[2]=-1.0;
      c[0]=-1.0;
      c[1]=-1.0;
      c[2]=-1.0;

      out[0] = coeff * (a[0] + b[0]*in[0] + b[1]*in[1] + b[2]*in[2]) * (a[1] + c[0]*in[0] + c[1]*in[1] + c[2]*in[2]);

      // case 1:
      coeff=2;
      a[0]=0.0;
      a[1]=-0.5;
      b[0]=1.0;
      b[1]=0.0;
      b[2]=0.0;
      c[0]=1.0;
      c[1]=0.0;
      c[2]=0.0;

      out[1] = coeff * (a[0] + b[0]*in[0] + b[1]*in[1] + b[2]*in[2]) * (a[1] + c[0]*in[0] + c[1]*in[1] + c[2]*in[2]);

      // case 2:
      coeff=2;
      a[0]=0.0;
      a[1]=-0.5;
      b[0]=0.0;
      b[1]=1.0;
      b[2]=0.0;
      c[0]=0.0;
      c[1]=1.0;
      c[2]=0.0;

      out[2] = coeff * (a[0] + b[0]*in[0] + b[1]*in[1] + b[2]*in[2]) * (a[1] + c[0]*in[0] + c[1]*in[1] + c[2]*in[2]);

      // case 3:
      coeff=2;
      a[0]=0.0;
      a[1]=-0.5;
      b[0]=0.0;
      b[1]=0.0;
      b[2]=1.0;
      c[0]=0.0;
      c[1]=0.0;
      c[2]=1.0;

      out[3] = coeff * (a[0] + b[0]*in[0] + b[1]*in[1] + b[2]*in[2]) * (a[1] + c[0]*in[0] + c[1]*in[1] + c[2]*in[2]);

      // case 4:
      coeff=4;
      a[0]=0.0;
      a[1]=1.0;
      b[0]=1.0;
      b[1]=0.0;
      b[2]=0.0;
      c[0]=-1.0;
      c[1]=-1.0;
      c[2]=-1.0;

      out[4] = coeff * (a[0] + b[0]*in[0] + b[1]*in[1] + b[2]*in[2]) * (a[1] + c[0]*in[0] + c[1]*in[1] + c[2]*in[2]);

      // case 5:
      coeff=4;
      a[0]=0.0;
      a[1]=0.0;
      b[0]=1.0;
      b[1]=0.0;
      b[2]=0.0;
      c[0]=0.0;
      c[1]=1.0;
      c[2]=0.0;

      out[5] = coeff * (a[0] + b[0]*in[0] + b[1]*in[1] + b[2]*in[2]) * (a[1] + c[0]*in[0] + c[1]*in[1] + c[2]*in[2]);

      // case 6:
      coeff=4;
      a[0]=0.0;
      a[1]=1.0;
      b[0]=0.0;
      b[1]=1.0;
      b[2]=0.0;
      c[0]=-1.0;
      c[1]=-1.0;
      c[2]=-1.0;

      out[6] = coeff * (a[0] + b[0]*in[0] + b[1]*in[1] + b[2]*in[2]) * (a[1] + c[0]*in[0] + c[1]*in[1] + c[2]*in[2]);

      // case 7:
      coeff=4;
      a[0]=0.0;
      a[1]=1.0;
      b[0]=0.0;
      b[1]=0.0;
      b[2]=1.0;
      c[0]=-1.0;
      c[1]=-1.0;
      c[2]=-1.0;

      out[7] = coeff * (a[0] + b[0]*in[0] + b[1]*in[1] + b[2]*in[2]) * (a[1] + c[0]*in[0] + c[1]*in[1] + c[2]*in[2]);

      // case 8:
      coeff=4;
      a[0]=0.0;
      a[1]=0.0;
      b[0]=1.0;
      b[1]=0.0;
      b[2]=0.0;
      c[0]=0.0;
      c[1]=0.0;
      c[2]=1.0;

      out[8] = coeff * (a[0] + b[0]*in[0] + b[1]*in[1] + b[2]*in[2]) * (a[1] + c[0]*in[0] + c[1]*in[1] + c[2]*in[2]);

      // case 9:
      coeff=4;
      a[0]=0.0;
      a[1]=0.0;
      b[0]=0.0;
      b[1]=1.0;
      b[2]=0.0;
      c[0]=0.0;
      c[1]=0.0;
      c[2]=1.0;

      out[9] = coeff * (a[0] + b[0]*in[0] + b[1]*in[1] + b[2]*in[2]) * (a[1] + c[0]*in[0] + c[1]*in[1] + c[2]*in[2]);

    }

    //! \brief Evaluate Jacobian of all shape functions
    inline void
    evaluateJacobian (const typename Traits::DomainType& in,         // position
                      std::vector<typename Traits::JacobianType>& out) const      // return value
    {
      out.resize(10);

      R aa[3][3], bb[3][3];
      // case 0:
      //x derivative
      aa[0][0]=-3.0;
      bb[0][0]=4.0;
      bb[1][0]=4.0;
      bb[2][0]=4.0;
      //y derivative
      aa[0][1]=-3.0;
      bb[0][1]=4.0;
      bb[1][1]=4.0;
      bb[2][1]=4.0;
      // z derivative
      aa[0][2]=-3.0;
      bb[0][2]=4.0;
      bb[1][2]=4.0;
      bb[2][2]=4.0;

      out[0][0][0] = aa[0][0] + bb[0][0]*in[0] + bb[1][0]*in[1] + bb[2][0]*in[2];
      out[0][0][1] = aa[0][1] + bb[0][1]*in[0] + bb[1][1]*in[1] + bb[2][1]*in[2];
      out[0][0][2] = aa[0][2] + bb[0][2]*in[0] + bb[1][2]*in[1] + bb[2][2]*in[2];

      //case 1:
      //x derivative
      aa[0][0]=-1.0;
      bb[0][0]=4.0;
      bb[1][0]=0.0;
      bb[2][0]=0.0;
      //y derivative
      aa[0][1]=0.0;
      bb[0][1]=0.0;
      bb[1][1]=0.0;
      bb[2][1]=0.0;
      // z derivative
      aa[0][2]=0.0;
      bb[0][2]=0.0;
      bb[1][2]=0.0;
      bb[2][2]=0.0;

      out[1][0][0] = aa[0][0] + bb[0][0]*in[0] + bb[1][0]*in[1] + bb[2][0]*in[2];
      out[1][0][1] = aa[0][1] + bb[0][1]*in[0] + bb[1][1]*in[1] + bb[2][1]*in[2];
      out[1][0][2] = aa[0][2] + bb[0][2]*in[0] + bb[1][2]*in[1] + bb[2][2]*in[2];

      // case 2:
      //x derivative
      aa[0][0]=0.0;
      bb[0][0]=0.0;
      bb[1][0]=0.0;
      bb[2][0]=0.0;
      //y derivative
      aa[0][1]=-1.0;
      bb[0][1]=0.0;
      bb[1][1]=4.0;
      bb[2][1]=0.0;
      // z derivative
      aa[0][2]=0.0;
      bb[0][2]=0.0;
      bb[1][2]=0.0;
      bb[2][2]=0.0;

      out[2][0][0] = aa[0][0] + bb[0][0]*in[0] + bb[1][0]*in[1] + bb[2][0]*in[2];
      out[2][0][1] = aa[0][1] + bb[0][1]*in[0] + bb[1][1]*in[1] + bb[2][1]*in[2];
      out[2][0][2] = aa[0][2] + bb[0][2]*in[0] + bb[1][2]*in[1] + bb[2][2]*in[2];

      // case 3:
      //x derivative
      aa[0][0]=0.0;
      bb[0][0]=0.0;
      bb[1][0]=0.0;
      bb[2][0]=0.0;
      //y derivative
      aa[0][1]=0.0;
      bb[0][1]=0.0;
      bb[1][1]=0.0;
      bb[2][1]=0.0;
      // z derivative
      aa[0][2]=-1.0;
      bb[0][2]=0.0;
      bb[1][2]=0.0;
      bb[2][2]=4.0;

      out[3][0][0] = aa[0][0] + bb[0][0]*in[0] + bb[1][0]*in[1] + bb[2][0]*in[2];
      out[3][0][1] = aa[0][1] + bb[0][1]*in[0] + bb[1][1]*in[1] + bb[2][1]*in[2];
      out[3][0][2] = aa[0][2] + bb[0][2]*in[0] + bb[1][2]*in[1] + bb[2][2]*in[2];

      // case 4:
      //x derivative
      aa[0][0]=4.0;
      bb[0][0]=-8.0;
      bb[1][0]=-4.0;
      bb[2][0]=-4.0;
      //y derivative
      aa[0][1]=0.0;
      bb[0][1]=-4.0;
      bb[1][1]=0.0;
      bb[2][1]=0.0;
      // z derivative
      aa[0][2]=0.0;
      bb[0][2]=-4.0;
      bb[1][2]=0.0;
      bb[2][2]=0.0;

      out[4][0][0] = aa[0][0] + bb[0][0]*in[0] + bb[1][0]*in[1] + bb[2][0]*in[2];
      out[4][0][1] = aa[0][1] + bb[0][1]*in[0] + bb[1][1]*in[1] + bb[2][1]*in[2];
      out[4][0][2] = aa[0][2] + bb[0][2]*in[0] + bb[1][2]*in[1] + bb[2][2]*in[2];

      // case 5:
      //x derivative
      aa[0][0]=0.0;
      bb[0][0]=0.0;
      bb[1][0]=4.0;
      bb[2][0]=0.0;
      //y derivative
      aa[0][1]=0.0;
      bb[0][1]=4.0;
      bb[1][1]=0.0;
      bb[2][1]=0.0;
      // z derivative
      aa[0][2]=0.0;
      bb[0][2]=0.0;
      bb[1][2]=0.0;
      bb[2][2]=0.0;

      out[5][0][0] = aa[0][0] + bb[0][0]*in[0] + bb[1][0]*in[1] + bb[2][0]*in[2];
      out[5][0][1] = aa[0][1] + bb[0][1]*in[0] + bb[1][1]*in[1] + bb[2][1]*in[2];
      out[5][0][2] = aa[0][2] + bb[0][2]*in[0] + bb[1][2]*in[1] + bb[2][2]*in[2];

      // case 6:
      //x derivative
      aa[0][0]=0.0;
      bb[0][0]=0.0;
      bb[1][0]=-4.0;
      bb[2][0]=0.0;
      //y derivative
      aa[0][1]=4.0;
      bb[0][1]=-4.0;
      bb[1][1]=-8.0;
      bb[2][1]=-4.0;
      // z derivative
      aa[0][2]=0.0;
      bb[0][2]=0.0;
      bb[1][2]=-4.0;
      bb[2][2]=0.0;

      out[6][0][0] = aa[0][0] + bb[0][0]*in[0] + bb[1][0]*in[1] + bb[2][0]*in[2];
      out[6][0][1] = aa[0][1] + bb[0][1]*in[0] + bb[1][1]*in[1] + bb[2][1]*in[2];
      out[6][0][2] = aa[0][2] + bb[0][2]*in[0] + bb[1][2]*in[1] + bb[2][2]*in[2];

      // case 7:
      //x derivative
      aa[0][0]=0.0;
      bb[0][0]=0.0;
      bb[1][0]=0.0;
      bb[2][0]=-4.0;
      //y derivative
      aa[0][1]=0.0;
      bb[0][1]=0.0;
      bb[1][1]=0.0;
      bb[2][1]=-4.0;
      // z derivative
      aa[0][2]=4.0;
      bb[0][2]=-4.0;
      bb[1][2]=-4.0;
      bb[2][2]=-8.0;

      out[7][0][0] = aa[0][0] + bb[0][0]*in[0] + bb[1][0]*in[1] + bb[2][0]*in[2];
      out[7][0][1] = aa[0][1] + bb[0][1]*in[0] + bb[1][1]*in[1] + bb[2][1]*in[2];
      out[7][0][2] = aa[0][2] + bb[0][2]*in[0] + bb[1][2]*in[1] + bb[2][2]*in[2];

      //case 8:
      //x derivative
      aa[0][0]=0.0;
      bb[0][0]=0.0;
      bb[1][0]=0.0;
      bb[2][0]=4.0;
      //y derivative
      aa[0][1]=0.0;
      bb[0][1]=0.0;
      bb[1][1]=0.0;
      bb[2][1]=0.0;
      // z derivative
      aa[0][2]=0.0;
      bb[0][2]=4.0;
      bb[1][2]=0.0;
      bb[2][2]=0.0;

      out[8][0][0] = aa[0][0] + bb[0][0]*in[0] + bb[1][0]*in[1] + bb[2][0]*in[2];
      out[8][0][1] = aa[0][1] + bb[0][1]*in[0] + bb[1][1]*in[1] + bb[2][1]*in[2];
      out[8][0][2] = aa[0][2] + bb[0][2]*in[0] + bb[1][2]*in[1] + bb[2][2]*in[2];

      // case 9:
      //x derivative
      aa[0][0]=0.0;
      bb[0][0]=0.0;
      bb[1][0]=0.0;
      bb[2][0]=0.0;
      //y derivative
      aa[0][1]=0.0;
      bb[0][1]=0.0;
      bb[1][1]=0.0;
      bb[2][1]=4.0;
      // z derivative
      aa[0][2]=0.0;
      bb[0][2]=0.0;
      bb[1][2]=4.0;
      bb[2][2]=0.0;

      out[9][0][0] = aa[0][0] + bb[0][0]*in[0] + bb[1][0]*in[1] + bb[2][0]*in[2];
      out[9][0][1] = aa[0][1] + bb[0][1]*in[0] + bb[1][1]*in[1] + bb[2][1]*in[2];
      out[9][0][2] = aa[0][2] + bb[0][2]*in[0] + bb[1][2]*in[1] + bb[2][2]*in[2];

    }

    /** \brief Evaluate partial derivatives of any order of all shape functions
     * \param order Order of the partial derivatives, in the classic multi-index notation
     * \param in Position where to evaluate the derivatives
     * \param[out] out Return value: the desired partial derivatives
     */
    void partial(const std::array<unsigned int,3>& order,
                 const typename Traits::DomainType& in,
                 std::vector<typename Traits::RangeType>& out) const
    {
      auto totalOrder = std::accumulate(order.begin(), order.end(), 0);
      if (totalOrder == 0) {
        evaluateFunction(in, out);
      } else if (totalOrder == 1) {
        out.resize(size());

        auto const direction = std::distance(order.begin(), std::find(order.begin(), order.end(), 1));
        switch (direction) {
          case 0:
            out[0] =-3.0 + 4.0*(in[0] + in[1] + in[2]);
            out[1] =-1.0 + 4.0*in[0];
            out[2] = 0.0;
            out[3] = 0.0;
            out[4] = 4.0 - 4.0*(2.0*in[0] + in[1] + in[2]);
            out[5] = 4.0*in[1];
            out[6] =-4.0*in[1];
            out[7] =-4.0*in[2];
            out[8] = 4.0*in[2];
            out[9] = 0.0;
            break;
          case 1:
            out[0] =-3.0 + 4.0*(in[0] + in[1] + in[2]);
            out[1] = 0.0;
            out[2] =-1.0 + 4.0*in[1];
            out[3] = 0.0;
            out[4] =-4.0*in[0];
            out[5] = 4.0*in[0];
            out[6] = 4.0 - 4.0*(in[0] + 2.0*in[1] + in[2]);
            out[7] =-4.0*in[2];
            out[8] = 0.0;
            out[9] = 4.0*in[2];
            break;
          case 2:
            out[0] =-3.0 + 4.0*(in[0] + in[1] + in[2]);
            out[1] = 0.0;
            out[2] = 0.0;
            out[3] =-1.0 + 4.0*in[2];
            out[4] =-4.0*in[0];
            out[5] = 0.0;
            out[6] =-4.0*in[1];
            out[7] = 4.0 - 4.0*(in[0] + in[1] + 2.0*in[2]);
            out[8] = 4.0*in[0];
            out[9] = 4.0*in[1];
            break;
          default:
            DUNE_THROW(RangeError, "Component out of range.");
        }
      } else {
        DUNE_THROW(NotImplemented, "Desired derivative order is not implemented");
      }
    }

    //! \brief Polynomial order of the shape functions
    unsigned int order () const
    {
      return 2;
    }
  };
}
#endif
