// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_Q2_LOCALBASIS_HH
#define DUNE_Q2_LOCALBASIS_HH

#include <dune/common/fmatrix.hh>

#include <dune/localfunctions/common/localbasis.hh>

namespace Dune
{
  /**@ingroup LocalBasisImplementation
     \brief Lagrange shape functions of order 2 on the reference cube

     Also known as \f$Q^2\f$.

     \tparam D Type to represent the field in the domain.
     \tparam R Type to represent the field in the range.
     \tparam dim Dimension of the reference cube

     \nosubgrouping
   */
  template<class D, class R, int dim>
  class Q2LocalBasis
  {
  public:
    typedef LocalBasisTraits<D,dim,Dune::FieldVector<D,dim>,R,1,Dune::FieldVector<R,1>,
        Dune::FieldMatrix<R,1,dim> > Traits;

    //! \brief number of shape functions
    unsigned int size () const
    {
      int size = 1;
      for (int i=0; i<dim; i++)
        size *= 3;
      return size;
    }

    //! \brief Evaluate all shape functions
    inline void evaluateFunction (const typename Traits::DomainType& in,
                                  std::vector<typename Traits::RangeType>& out) const
    {
      out.resize(size());

      // Evaluate the Lagrange functions
      array<array<R,3>, dim> X;

      for (size_t i=0; i<dim; i++) {
        X[i][0] =  R(2)*in[i]*in[i] - R(3)*in[i]+R(1);
        X[i][1] = -R(4)*in[i]*in[i] + R(4)*in[i];
        X[i][2] =  R(2)*in[i]*in[i] -   in[i];
      }

      // legacy special case
      if (dim==3) {

        R x=in[0], y=in[1], z=in[2];
        R X0=R(2)*x*x-R(3)*x+R(1), X1=-R(4)*x*x+R(4)*x, X2=R(2)*x*x-x;
        R Y0=R(2)*y*y-R(3)*y+R(1), Y1=-R(4)*y*y+R(4)*y, Y2=R(2)*y*y-y;
        R Z0=R(2)*z*z-R(3)*z+R(1), Z1=-R(4)*z*z+R(4)*z, Z2=R(2)*z*z-z;

        // numbering: first in descending codim order,
        // second according to the reference element numbering
        out[0]  = X0*Y0*Z0;
        out[14] = X1*Y0*Z0;
        out[1]  = X2*Y0*Z0;
        out[12] = X0*Y1*Z0;
        out[24] = X1*Y1*Z0;
        out[13] = X2*Y1*Z0;
        out[2]  = X0*Y2*Z0;
        out[15] = X1*Y2*Z0;
        out[3]  = X2*Y2*Z0;

        out[8]  = X0*Y0*Z1;
        out[22] = X1*Y0*Z1;
        out[9]  = X2*Y0*Z1;
        out[20] = X0*Y1*Z1;
        out[26] = X1*Y1*Z1;
        out[21] = X2*Y1*Z1;
        out[10] = X0*Y2*Z1;
        out[23] = X1*Y2*Z1;
        out[11] = X2*Y2*Z1;

        out[4]  = X0*Y0*Z2;
        out[18] = X1*Y0*Z2;
        out[5]  = X2*Y0*Z2;
        out[16] = X0*Y1*Z2;
        out[25] = X1*Y1*Z2;
        out[17] = X2*Y1*Z2;
        out[6]  = X0*Y2*Z2;
        out[19] = X1*Y2*Z2;
        out[7]  = X2*Y2*Z2;

        return;
      }

      for (size_t i=0; i<out.size(); i++) {

        out[i] = 1;

        // Construct the i-th Lagrange point
        size_t ternary = i;
        for (int j=0; j<dim; j++) {

          int digit = ternary%3;
          ternary /= 3;

          // Multiply the 1d Lagrange shape functions together
          out[i] *= X[j][digit];

        }

      }

    }

    //! \brief Evaluate Jacobian of all shape functions
    inline void
    evaluateJacobian (const typename Traits::DomainType& in,       // position
                      std::vector<typename Traits::JacobianType>& out) const // return value
    {
      out.resize(size());

      switch (dim) {

      case 1 : {

        out[0][0][0] =  R(4)*in[0] - R(3);
        out[1][0][0] = -R(8)*in[0] + R(4);
        out[2][0][0] =  R(4)*in[0] - R(1);
        break;
      }

      case 2 : {

        R x=in[0], y=in[1];
        R X0=R(2)*x*x-R(3)*x+R(1), X1=-R(4)*x*x+R(4)*x, X2=R(2)*x*x-x;
        R Y0=R(2)*y*y-R(3)*y+R(1), Y1=-R(4)*y*y+R(4)*y, Y2=R(2)*y*y-y;
        R DX0=R(4)*x-R(3), DX1=-R(8)*x+R(4), DX2=R(4)*x-R(1);
        R DY0=R(4)*y-R(3), DY1=-R(8)*y+R(4), DY2=R(4)*y-R(1);

        out[6][0][0] = DX0*Y2; out[7][0][0] = DX1*Y2; out[8][0][0] = DX2*Y2;
        out[6][0][1] = X0*DY2; out[7][0][1] = X1*DY2; out[8][0][1] = X2*DY2;

        out[3][0][0] = DX0*Y1; out[4][0][0] = DX1*Y1; out[5][0][0] = DX2*Y1;
        out[3][0][1] = X0*DY1; out[4][0][1] = X1*DY1; out[5][0][1] = X2*DY1;

        out[0][0][0] = DX0*Y0; out[1][0][0] = DX1*Y0; out[2][0][0] = DX2*Y0;
        out[0][0][1] = X0*DY0; out[1][0][1] = X1*DY0; out[2][0][1] = X2*DY0;
        break;
      }

      case 3 : {
        R x=in[0], y=in[1], z=in[2];
        R X0=R(2)*x*x-R(3)*x+R(1), X1=-R(4)*x*x+R(4)*x, X2=R(2)*x*x-x;
        R Y0=R(2)*y*y-R(3)*y+R(1), Y1=-R(4)*y*y+R(4)*y, Y2=R(2)*y*y-y;
        R Z0=R(2)*z*z-R(3)*z+R(1), Z1=-R(4)*z*z+R(4)*z, Z2=R(2)*z*z-z;
        R DX0=R(4)*x-R(3), DX1=-R(8)*x+R(4), DX2=R(4)*x-R(1);
        R DY0=R(4)*y-R(3), DY1=-R(8)*y+R(4), DY2=R(4)*y-R(1);
        R DZ0=R(4)*z-R(3), DZ1=-R(8)*z+R(4), DZ2=R(4)*z-R(1);

        out[0][0][0]  = DX0*Y0*Z0; out[0][0][1]  = X0*DY0*Z0; out[0][0][2]  = X0*Y0*DZ0;
        out[14][0][0] = DX1*Y0*Z0; out[14][0][1] = X1*DY0*Z0; out[14][0][2] = X1*Y0*DZ0;
        out[1][0][0]  = DX2*Y0*Z0; out[1][0][1]  = X2*DY0*Z0; out[1][0][2]  = X2*Y0*DZ0;
        out[12][0][0] = DX0*Y1*Z0; out[12][0][1] = X0*DY1*Z0; out[12][0][2] = X0*Y1*DZ0;
        out[24][0][0] = DX1*Y1*Z0; out[24][0][1] = X1*DY1*Z0; out[24][0][2] = X1*Y1*DZ0;
        out[13][0][0] = DX2*Y1*Z0; out[13][0][1] = X2*DY1*Z0; out[13][0][2] = X2*Y1*DZ0;
        out[2][0][0]  = DX0*Y2*Z0; out[2][0][1]  = X0*DY2*Z0; out[2][0][2]  = X0*Y2*DZ0;
        out[15][0][0] = DX1*Y2*Z0; out[15][0][1] = X1*DY2*Z0; out[15][0][2] = X1*Y2*DZ0;
        out[3][0][0]  = DX2*Y2*Z0; out[3][0][1]  = X2*DY2*Z0; out[3][0][2]  = X2*Y2*DZ0;

        out[8][0][0]  = DX0*Y0*Z1; out[8][0][1]  = X0*DY0*Z1; out[8][0][2]  = X0*Y0*DZ1;
        out[22][0][0] = DX1*Y0*Z1; out[22][0][1] = X1*DY0*Z1; out[22][0][2] = X1*Y0*DZ1;
        out[9][0][0]  = DX2*Y0*Z1; out[9][0][1]  = X2*DY0*Z1; out[9][0][2]  = X2*Y0*DZ1;
        out[20][0][0] = DX0*Y1*Z1; out[20][0][1] = X0*DY1*Z1; out[20][0][2] = X0*Y1*DZ1;
        out[26][0][0] = DX1*Y1*Z1; out[26][0][1] = X1*DY1*Z1; out[26][0][2] = X1*Y1*DZ1;
        out[21][0][0] = DX2*Y1*Z1; out[21][0][1] = X2*DY1*Z1; out[21][0][2] = X2*Y1*DZ1;
        out[10][0][0] = DX0*Y2*Z1; out[10][0][1] = X0*DY2*Z1; out[10][0][2] = X0*Y2*DZ1;
        out[23][0][0] = DX1*Y2*Z1; out[23][0][1] = X1*DY2*Z1; out[23][0][2] = X1*Y2*DZ1;
        out[11][0][0] = DX2*Y2*Z1; out[11][0][1] = X2*DY2*Z1; out[11][0][2] = X2*Y2*DZ1;

        out[4][0][0]  = DX0*Y0*Z2; out[4][0][1]  = X0*DY0*Z2; out[4][0][2]  = X0*Y0*DZ2;
        out[18][0][0] = DX1*Y0*Z2; out[18][0][1] = X1*DY0*Z2; out[18][0][2] = X1*Y0*DZ2;
        out[5][0][0]  = DX2*Y0*Z2; out[5][0][1]  = X2*DY0*Z2; out[5][0][2]  = X2*Y0*DZ2;
        out[16][0][0] = DX0*Y1*Z2; out[16][0][1] = X0*DY1*Z2; out[16][0][2] = X0*Y1*DZ2;
        out[25][0][0] = DX1*Y1*Z2; out[25][0][1] = X1*DY1*Z2; out[25][0][2] = X1*Y1*DZ2;
        out[17][0][0] = DX2*Y1*Z2; out[17][0][1] = X2*DY1*Z2; out[17][0][2] = X2*Y1*DZ2;
        out[6][0][0]  = DX0*Y2*Z2; out[6][0][1]  = X0*DY2*Z2; out[6][0][2]  = X0*Y2*DZ2;
        out[19][0][0] = DX1*Y2*Z2; out[19][0][1] = X1*DY2*Z2; out[19][0][2] = X1*Y2*DZ2;
        out[7][0][0]  = DX2*Y2*Z2; out[7][0][1]  = X2*DY2*Z2; out[7][0][2]  = X2*Y2*DZ2;

        break;
      }
      default :
        DUNE_THROW(NotImplemented, "Q2LocalBasis for dim==" << dim);

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
