// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_Q22DLOCALBASIS_HH
#define DUNE_Q22DLOCALBASIS_HH

#include "../common/localbasis.hh"

namespace Dune
{
  /**@ingroup LocalBasisImplementation
     \brief Lagrange shape functions of order 2 on the reference quadrilateral.

     Also known as \f$Q^2\f$.

     \tparam D Type to represent the field in the domain.
     \tparam R Type to represent the field in the range.

     \nosubgrouping
   */
  template<class D, class R>
  class Q22DLocalBasis
  {
  public:
    typedef C1LocalBasisTraits<D,2,Dune::FieldVector<D,2>,R,1,Dune::FieldVector<R,1>,
        Dune::FieldVector<Dune::FieldVector<R,2>,1> > Traits;

    //! \brief number of shape functions
    unsigned int size () const
    {
      return 9;
    }

    //! \brief Evaluate all shape functions
    inline void evaluateFunction (const typename Traits::DomainType& in,
                                  std::vector<typename Traits::RangeType>& out) const
    {
      out.resize(9);

      R x=in[0], y=in[1];
      R X0=2*x*x-3*x+1, X1=-4*x*x+4*x, X2=2*x*x-x;
      R Y0=2*y*y-3*y+1, Y1=-4*y*y+4*y, Y2=2*y*y-y;

      out[2] = X0*Y2;
      out[7] = X1*Y2;
      out[3] = X2*Y2;

      out[4] = X0*Y1;
      out[8] = X1*Y1;
      out[5] = X2*Y1;

      out[0] = X0*Y0;
      out[6] = X1*Y0;
      out[1] = X2*Y0;
    }

    //! \brief Evaluate Jacobian of all shape functions
    inline void
    evaluateJacobian (const typename Traits::DomainType& in,       // position
                      std::vector<typename Traits::JacobianType>& out) const        // return value
    {
      out.resize(9);

      R x=in[0], y=in[1];
      R X0=2*x*x-3*x+1, X1=4*x*x-4*x, X2=2*x*x-x;
      R Y0=2*y*y-3*y+1, Y1=4*y*y-4*y, Y2=2*y*y-y;
      R DX0=4*x-3, DX1=-8*x+4, DX2=4*x-1;
      R DY0=4*y-3, DY1=-8*y+4, DY2=4*y-1;

      out[2][0][0] = DX0*Y2; out[7][0][0] = DX1*Y2; out[3][0][0] = DX2*Y2;
      out[2][0][1] = X0*DY2; out[7][0][1] = X1*DY2; out[3][0][1] = X2*DY2;

      out[4][0][0] = DX0*Y1; out[8][0][0] = DX1*Y1; out[5][0][0] = DX2*Y1;
      out[4][0][1] = X0*DY1; out[8][0][1] = X1*DY1; out[5][0][1] = X2*DY1;

      out[0][0][0] = DX0*Y0; out[6][0][0] = DX1*Y0; out[1][0][0] = DX2*Y0;
      out[0][0][1] = X0*DY0; out[6][0][1] = X1*DY0; out[1][0][1] = X2*DY0;
    }

    //! \brief Polynomial order of the shape functions
    unsigned int order () const
    {
      return 2;
    }
  };
}
#endif
