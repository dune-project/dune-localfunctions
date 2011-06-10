// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_Q23DLOCALBASIS_HH
#define DUNE_Q23DLOCALBASIS_HH

#include <dune/common/fmatrix.hh>

#include <dune/localfunctions/common/localbasis.hh>

namespace Dune
{
  /**@ingroup LocalBasisImplementation
     \brief Lagrange shape functions of order 2 on the reference hexahedron.

     Also known as \f$Q^2\f$.

     \tparam D Type to represent the field in the domain.
     \tparam R Type to represent the field in the range.

     \nosubgrouping
   */
  template<class D, class R>
  class Q23DLocalBasis
  {
  public:
    typedef LocalBasisTraits<D,3,Dune::FieldVector<D,3>,R,1,Dune::FieldVector<R,1>,
        Dune::FieldMatrix<R,1,3> > Traits;

    //! \brief number of shape functions
    unsigned int size () const
    {
      return 27;
    }

    //! \brief Evaluate all shape functions
    inline void evaluateFunction (const typename Traits::DomainType& in,
                                  std::vector<typename Traits::RangeType>& out) const
    {
      out.resize(27);

      R x=in[0], y=in[1], z=in[2];
      R X0=2*x*x-3*x+1, X1=-4*x*x+4*x, X2=2*x*x-x;
      R Y0=2*y*y-3*y+1, Y1=-4*y*y+4*y, Y2=2*y*y-y;
      R Z0=2*z*z-3*z+1, Z1=-4*z*z+4*z, Z2=2*z*z-z;

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
    }

    //! \brief Evaluate Jacobian of all shape functions
    inline void
    evaluateJacobian (const typename Traits::DomainType& in,       // position
                      std::vector<typename Traits::JacobianType>& out) const // return value
    {
      out.resize(27);

      R x=in[0], y=in[1], z=in[2];
      R X0=2*x*x-3*x+1, X1=-4*x*x+4*x, X2=2*x*x-x;
      R Y0=2*y*y-3*y+1, Y1=-4*y*y+4*y, Y2=2*y*y-y;
      R Z0=2*z*z-3*z+1, Z1=-4*z*z+4*z, Z2=2*z*z-z;
      R DX0=4*x-3, DX1=-8*x+4, DX2=4*x-1;
      R DY0=4*y-3, DY1=-8*y+4, DY2=4*y-1;
      R DZ0=4*z-3, DZ1=-8*z+4, DZ2=4*z-1;

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
    }

    //! \brief Polynomial order of the shape functions
    unsigned int order () const
    {
      return 2;
    }
  };
}
#endif
