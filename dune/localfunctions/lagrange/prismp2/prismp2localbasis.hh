// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PRISM_P2_LOCALBASIS_HH
#define DUNE_PRISM_P2_LOCALBASIS_HH

#include <numeric>

#include <dune/common/fmatrix.hh>

#include <dune/localfunctions/common/localbasis.hh>

namespace Dune
{
  /**@ingroup LocalBasisImplementation
         \brief Quadratic Lagrange shape functions on the prism.

         Defines the quadratic shape functions on prism.

         \tparam D Type to represent the field in the domain.
         \tparam R Type to represent the field in the range.

         \nosubgrouping
   */
  template<class D, class R>
  class PrismP2LocalBasis
  {
  public:
    //! \brief export type traits for function signature
    typedef LocalBasisTraits<D,3,Dune::FieldVector<D,3>,R,1,Dune::FieldVector<R,1>,
        Dune::FieldMatrix<R,1,3> > Traits;

    //! \brief number of shape functions
    unsigned int size () const
    {
      return 18;
    }

    //! \brief Evaluate all shape functions
    inline void evaluateFunction (const typename Traits::DomainType& in,
                                  std::vector<typename Traits::RangeType>& out) const
    {
      out.resize(18);


      int coeff;
      R a[2], b[2], c[2], a1d, b1d, c1d;


      // lower triangle:
      coeff= 2;
      a[0] = 1;
      a[1] = 0.5;
      b[0] = -1;
      b[1] = -1;
      c[0] = -1;
      c[1] = -1;
      a1d  = 1;
      b1d  = -3;
      c1d  = 2;
      out[0] = coeff * (a[0] + b[0]*in[0] + b[1]*in[1]) * (a[1] + c[0]*in[0] + c[1]*in[1])*(a1d + in[2]*b1d + in[2]*in[2]*c1d);


      coeff= 2;
      a[0] = 0;
      a[1] = -0.5;
      b[0] = 1;
      b[1] = 0;
      c[0] = 1;
      c[1] = 0;
      a1d  = 1;
      b1d  = -3;
      c1d  = 2;
      out[1] = coeff * (a[0] + b[0]*in[0] + b[1]*in[1]) * (a[1] + c[0]*in[0] + c[1]*in[1])*(a1d + in[2]*b1d + in[2]*in[2]*c1d);


      coeff= 2;
      a[0] = 0;
      a[1] = -0.5;
      b[0] = 0;
      b[1] = 1;
      c[0] = 0;
      c[1] = 1;
      a1d  = 1;
      b1d  = -3;
      c1d  = 2;
      out[2] = coeff * (a[0] + b[0]*in[0] + b[1]*in[1]) * (a[1] + c[0]*in[0] + c[1]*in[1])*(a1d + in[2]*b1d + in[2]*in[2]*c1d);

      //upper triangle
      coeff= 2;
      a[0] = 1;
      a[1] = 0.5;
      b[0] = -1;
      b[1] = -1;
      c[0] = -1;
      c[1] = -1;
      a1d  = 0;
      b1d  = -1;
      c1d  = 2;
      out[3] = coeff * (a[0] + b[0]*in[0] + b[1]*in[1]) * (a[1] + c[0]*in[0] + c[1]*in[1])*(a1d + in[2]*b1d + in[2]*in[2]*c1d);


      coeff= 2;
      a[0] = 0;
      a[1] = -0.5;
      b[0] = 1;
      b[1] = 0;
      c[0] = 1;
      c[1] = 0;
      a1d  = 0;
      b1d  = -1;
      c1d  = 2;
      out[4] = coeff * (a[0] + b[0]*in[0] + b[1]*in[1]) * (a[1] + c[0]*in[0] + c[1]*in[1])*(a1d + in[2]*b1d + in[2]*in[2]*c1d);


      coeff= 2;
      a[0] = 0;
      a[1] = -0.5;
      b[0] = 0;
      b[1] = 1;
      c[0] = 0;
      c[1] = 1;
      a1d  = 0;
      b1d  = -1;
      c1d  = 2;
      out[5] = coeff * (a[0] + b[0]*in[0] + b[1]*in[1]) * (a[1] + c[0]*in[0] + c[1]*in[1])*(a1d + in[2]*b1d + in[2]*in[2]*c1d);

      // vertical edges
      coeff= 2;
      a[0] = 1;
      a[1] = 0.5;
      b[0] = -1;
      b[1] = -1;
      c[0] = -1;
      c[1] = -1;
      a1d  = 0;
      b1d  = 4;
      c1d  = -4;
      out[6] = coeff * (a[0] + b[0]*in[0] + b[1]*in[1]) * (a[1] + c[0]*in[0] + c[1]*in[1])*(a1d + in[2]*b1d + in[2]*in[2]*c1d);


      coeff= 2;
      a[0] = 0;
      a[1] = -0.5;
      b[0] = 1;
      b[1] = 0;
      c[0] = 1;
      c[1] = 0;
      a1d  = 0;
      b1d  = 4;
      c1d  = -4;
      out[7] = coeff * (a[0] + b[0]*in[0] + b[1]*in[1]) * (a[1] + c[0]*in[0] + c[1]*in[1])*(a1d + in[2]*b1d + in[2]*in[2]*c1d);


      coeff= 2;
      a[0] = 0;
      a[1] = -0.5;
      b[0] = 0;
      b[1] = 1;
      c[0] = 0;
      c[1] = 1;
      a1d  = 0;
      b1d  = 4;
      c1d  = -4;
      out[8] = coeff * (a[0] + b[0]*in[0] + b[1]*in[1]) * (a[1] + c[0]*in[0] + c[1]*in[1])*(a1d + in[2]*b1d + in[2]*in[2]*c1d);

      // lower triangle edges
      coeff= 4;
      a[0] = 0;
      a[1] = 1;
      b[0] = 1;
      b[1] = 0;
      c[0] = -1;
      c[1] = -1;
      a1d  = 1;
      b1d  = -3;
      c1d  = 2;
      out[9] = coeff * (a[0] + b[0]*in[0] + b[1]*in[1]) * (a[1] + c[0]*in[0] + c[1]*in[1])*(a1d + in[2]*b1d + in[2]*in[2]*c1d);


      coeff= 4;
      a[0] = 0;
      a[1] = 1;
      b[0] = 0;
      b[1] = 1;
      c[0] = -1;
      c[1] = -1;
      a1d  = 1;
      b1d  = -3;
      c1d  = 2;
      out[10] = coeff * (a[0] + b[0]*in[0] + b[1]*in[1]) * (a[1] + c[0]*in[0] + c[1]*in[1])*(a1d + in[2]*b1d + in[2]*in[2]*c1d);


      coeff= 4;
      a[0] = 0;
      a[1] = 0;
      b[0] = 1;
      b[1] = 0;
      c[0] = 0;
      c[1] = 1;
      a1d  = 1;
      b1d  = -3;
      c1d  = 2;
      out[11] = coeff * (a[0] + b[0]*in[0] + b[1]*in[1]) * (a[1] + c[0]*in[0] + c[1]*in[1])*(a1d + in[2]*b1d + in[2]*in[2]*c1d);

      // upper triangle edges
      coeff= 4;
      a[0] = 0;
      a[1] = 1;
      b[0] = 1;
      b[1] = 0;
      c[0] = -1;
      c[1] = -1;
      a1d  = 0;
      b1d  = -1;
      c1d  = 2;
      out[12] = coeff * (a[0] + b[0]*in[0] + b[1]*in[1]) * (a[1] + c[0]*in[0] + c[1]*in[1])*(a1d + in[2]*b1d + in[2]*in[2]*c1d);


      coeff= 4;
      a[0] = 0;
      a[1] = 1;
      b[0] = 0;
      b[1] = 1;
      c[0] = -1;
      c[1] = -1;
      a1d  = 0;
      b1d  = -1;
      c1d  = 2;
      out[13] = coeff * (a[0] + b[0]*in[0] + b[1]*in[1]) * (a[1] + c[0]*in[0] + c[1]*in[1])*(a1d + in[2]*b1d + in[2]*in[2]*c1d);


      coeff= 4;
      a[0] = 0;
      a[1] = 0;
      b[0] = 1;
      b[1] = 0;
      c[0] = 0;
      c[1] = 1;
      a1d  = 0;
      b1d  = -1;
      c1d  = 2;
      out[14] = coeff * (a[0] + b[0]*in[0] + b[1]*in[1]) * (a[1] + c[0]*in[0] + c[1]*in[1])*(a1d + in[2]*b1d + in[2]*in[2]*c1d);

      // quadrilateral sides
      coeff= 4;
      a[0] = 0;
      a[1] = 1;
      b[0] = 1;
      b[1] = 0;
      c[0] = -1;
      c[1] = -1;
      a1d  = 0;
      b1d  = 4;
      c1d  = -4;
      out[15] = coeff * (a[0] + b[0]*in[0] + b[1]*in[1]) * (a[1] + c[0]*in[0] + c[1]*in[1])*(a1d + in[2]*b1d + in[2]*in[2]*c1d);


      coeff= 4;
      a[0] = 0;
      a[1] = 1;
      b[0] = 0;
      b[1] = 1;
      c[0] = -1;
      c[1] = -1;
      a1d  = 0;
      b1d  = 4;
      c1d  = -4;
      out[16] = coeff * (a[0] + b[0]*in[0] + b[1]*in[1]) * (a[1] + c[0]*in[0] + c[1]*in[1])*(a1d + in[2]*b1d + in[2]*in[2]*c1d);


      coeff= 4;
      a[0] = 0;
      a[1] = 0;
      b[0] = 1;
      b[1] = 0;
      c[0] = 0;
      c[1] = 1;
      a1d  = 0;
      b1d  = 4;
      c1d  = -4;
      out[17] = coeff * (a[0] + b[0]*in[0] + b[1]*in[1]) * (a[1] + c[0]*in[0] + c[1]*in[1])*(a1d + in[2]*b1d + in[2]*in[2]*c1d);
    }

    //! \brief Evaluate Jacobian of all shape functions
    inline void
    evaluateJacobian (const typename Traits::DomainType& in,         // position
                      std::vector<typename Traits::JacobianType>& out) const      // return value
    {
      out.resize(18);



      int coeff;
      R a[2], b[2], c[2], aa[2], bb[2][2], a1d, b1d, c1d;


      // lower triangle:
      coeff= 2;
      a[0] = 1;
      a[1] = 0.5;
      b[0] = -1;
      b[1] = -1;
      c[0] = -1;
      c[1] = -1;
      a1d  = 1;
      b1d  = -3;
      c1d  = 2;
      aa[0] = -3;
      aa[1] = -3;
      bb[0][0] = 4;
      bb[0][1] = 4;
      bb[1][0] = 4;
      bb[1][1] = 4;
      out[0][0][0] = (aa[0] + bb[0][0]*in[0] + bb[1][0]*in[1]) * (a1d + in[2]*b1d + in[2]*in[2]*c1d);
      out[0][0][1] = (aa[1] + bb[0][1]*in[0] + bb[1][1]*in[1]) * (a1d + in[2]*b1d + in[2]*in[2]*c1d);
      out[0][0][2] = coeff * (a[0] + b[0]*in[0] + b[1]*in[1]) * (a[1] + c[0]*in[0] + c[1]*in[1]) * (b1d + 2*c1d*in[2]);



      coeff= 2;
      a[0] = 0;
      a[1] = -0.5;
      b[0] = 1;
      b[1] = 0;
      c[0] = 1;
      c[1] = 0;
      a1d  = 1;
      b1d  = -3;
      c1d  = 2;
      aa[0] = -1;
      aa[1] = 0;
      bb[0][0] = 4;
      bb[0][1] = 0;
      bb[1][0] = 0;
      bb[1][1] = 0;
      out[1][0][0] = (aa[0] + bb[0][0]*in[0] + bb[1][0]*in[1]) * (a1d + in[2]*b1d + in[2]*in[2]*c1d);
      out[1][0][1] = (aa[1] + bb[0][1]*in[0] + bb[1][1]*in[1]) * (a1d + in[2]*b1d + in[2]*in[2]*c1d);
      out[1][0][2] = coeff * (a[0] + b[0]*in[0] + b[1]*in[1]) * (a[1] + c[0]*in[0] + c[1]*in[1]) * (b1d + 2*c1d*in[2]);


      coeff= 2;
      a[0] = 0;
      a[1] = -0.5;
      b[0] = 0;
      b[1] = 1;
      c[0] = 0;
      c[1] = 1;
      a1d  = 1;
      b1d  = -3;
      c1d  = 2;
      aa[0] = 0;
      aa[1] = -1;
      bb[0][0] = 0;
      bb[0][1] = 0;
      bb[1][0] = 0;
      bb[1][1] = 4;
      out[2][0][0] = (aa[0] + bb[0][0]*in[0] + bb[1][0]*in[1]) * (a1d + in[2]*b1d + in[2]*in[2]*c1d);
      out[2][0][1] = (aa[1] + bb[0][1]*in[0] + bb[1][1]*in[1]) * (a1d + in[2]*b1d + in[2]*in[2]*c1d);
      out[2][0][2] = coeff * (a[0] + b[0]*in[0] + b[1]*in[1]) * (a[1] + c[0]*in[0] + c[1]*in[1]) * (b1d + 2*c1d*in[2]);


      //upper triangle
      coeff= 2;
      a[0] = 1;
      a[1] = 0.5;
      b[0] = -1;
      b[1] = -1;
      c[0] = -1;
      c[1] = -1;
      a1d  = 0;
      b1d  = -1;
      c1d  = 2;
      aa[0] = -3;
      aa[1] = -3;
      bb[0][0] = 4;
      bb[0][1] = 4;
      bb[1][0] = 4;
      bb[1][1] = 4;
      out[3][0][0] = (aa[0] + bb[0][0]*in[0] + bb[1][0]*in[1]) * (a1d + in[2]*b1d + in[2]*in[2]*c1d);
      out[3][0][1] = (aa[1] + bb[0][1]*in[0] + bb[1][1]*in[1]) * (a1d + in[2]*b1d + in[2]*in[2]*c1d);
      out[3][0][2] = coeff * (a[0] + b[0]*in[0] + b[1]*in[1]) * (a[1] + c[0]*in[0] + c[1]*in[1]) * (b1d + 2*c1d*in[2]);



      coeff= 2;
      a[0] = 0;
      a[1] = -0.5;
      b[0] = 1;
      b[1] = 0;
      c[0] = 1;
      c[1] = 0;
      a1d  = 0;
      b1d  = -1;
      c1d  = 2;
      aa[0] = -1;
      aa[1] = 0;
      bb[0][0] = 4;
      bb[0][1] = 0;
      bb[1][0] = 0;
      bb[1][1] = 0;
      out[4][0][0] = (aa[0] + bb[0][0]*in[0] + bb[1][0]*in[1]) * (a1d + in[2]*b1d + in[2]*in[2]*c1d);
      out[4][0][1] = (aa[1] + bb[0][1]*in[0] + bb[1][1]*in[1]) * (a1d + in[2]*b1d + in[2]*in[2]*c1d);
      out[4][0][2] = coeff * (a[0] + b[0]*in[0] + b[1]*in[1]) * (a[1] + c[0]*in[0] + c[1]*in[1]) * (b1d + 2*c1d*in[2]);



      coeff= 2;
      a[0] = 0;
      a[1] = -0.5;
      b[0] = 0;
      b[1] = 1;
      c[0] = 0;
      c[1] = 1;
      a1d  = 0;
      b1d  = -1;
      c1d  = 2;
      aa[0] = 0;
      aa[1] = -1;
      bb[0][0] = 0;
      bb[0][1] = 0;
      bb[1][0] = 0;
      bb[1][1] = 4;
      out[5][0][0] = (aa[0] + bb[0][0]*in[0] + bb[1][0]*in[1]) * (a1d + in[2]*b1d + in[2]*in[2]*c1d);
      out[5][0][1] = (aa[1] + bb[0][1]*in[0] + bb[1][1]*in[1]) * (a1d + in[2]*b1d + in[2]*in[2]*c1d);
      out[5][0][2] = coeff * (a[0] + b[0]*in[0] + b[1]*in[1]) * (a[1] + c[0]*in[0] + c[1]*in[1]) * (b1d + 2*c1d*in[2]);



      // vertical edges
      coeff= 2;
      a[0] = 1;
      a[1] = 0.5;
      b[0] = -1;
      b[1] = -1;
      c[0] = -1;
      c[1] = -1;
      a1d  = 0;
      b1d  = 4;
      c1d  = -4;
      aa[0] = -3;
      aa[1] = -3;
      bb[0][0] = 4;
      bb[0][1] = 4;
      bb[1][0] = 4;
      bb[1][1] = 4;
      out[6][0][0] = (aa[0] + bb[0][0]*in[0] + bb[1][0]*in[1]) * (a1d + in[2]*b1d + in[2]*in[2]*c1d);
      out[6][0][1] = (aa[1] + bb[0][1]*in[0] + bb[1][1]*in[1]) * (a1d + in[2]*b1d + in[2]*in[2]*c1d);
      out[6][0][2] = coeff * (a[0] + b[0]*in[0] + b[1]*in[1]) * (a[1] + c[0]*in[0] + c[1]*in[1]) * (b1d + 2*c1d*in[2]);



      coeff= 2;
      a[0] = 0;
      a[1] = -0.5;
      b[0] = 1;
      b[1] = 0;
      c[0] = 1;
      c[1] = 0;
      a1d  = 0;
      b1d  = 4;
      c1d  = -4;
      aa[0] = -1;
      aa[1] = 0;
      bb[0][0] = 4;
      bb[0][1] = 0;
      bb[1][0] = 0;
      bb[1][1] = 0;
      out[7][0][0] = (aa[0] + bb[0][0]*in[0] + bb[1][0]*in[1]) * (a1d + in[2]*b1d + in[2]*in[2]*c1d);
      out[7][0][1] = (aa[1] + bb[0][1]*in[0] + bb[1][1]*in[1]) * (a1d + in[2]*b1d + in[2]*in[2]*c1d);
      out[7][0][2] = coeff * (a[0] + b[0]*in[0] + b[1]*in[1]) * (a[1] + c[0]*in[0] + c[1]*in[1]) * (b1d + 2*c1d*in[2]);



      coeff= 2;
      a[0] = 0;
      a[1] = -0.5;
      b[0] = 0;
      b[1] = 1;
      c[0] = 0;
      c[1] = 1;
      a1d  = 0;
      b1d  = 4;
      c1d  = -4;
      aa[0] = 0;
      aa[1] = -1;
      bb[0][0] = 0;
      bb[0][1] = 0;
      bb[1][0] = 0;
      bb[1][1] = 4;
      out[8][0][0] = (aa[0] + bb[0][0]*in[0] + bb[1][0]*in[1]) * (a1d + in[2]*b1d + in[2]*in[2]*c1d);
      out[8][0][1] = (aa[1] + bb[0][1]*in[0] + bb[1][1]*in[1]) * (a1d + in[2]*b1d + in[2]*in[2]*c1d);
      out[8][0][2] = coeff * (a[0] + b[0]*in[0] + b[1]*in[1]) * (a[1] + c[0]*in[0] + c[1]*in[1]) * (b1d + 2*c1d*in[2]);




      // lower triangle edges
      coeff= 4;
      a[0] = 0;
      a[1] = 1;
      b[0] = 1;
      b[1] = 0;
      c[0] = -1;
      c[1] = -1;
      a1d  = 1;
      b1d  = -3;
      c1d  = 2;
      aa[0] = 4;
      aa[1] = 0;
      bb[0][0] = -8;
      bb[0][1] = -4;
      bb[1][0] = -4;
      bb[1][1] = 0;
      out[9][0][0] = (aa[0] + bb[0][0]*in[0] + bb[1][0]*in[1]) * (a1d + in[2]*b1d + in[2]*in[2]*c1d);
      out[9][0][1] = (aa[1] + bb[0][1]*in[0] + bb[1][1]*in[1]) * (a1d + in[2]*b1d + in[2]*in[2]*c1d);
      out[9][0][2] = coeff * (a[0] + b[0]*in[0] + b[1]*in[1]) * (a[1] + c[0]*in[0] + c[1]*in[1]) * (b1d + 2*c1d*in[2]);




      coeff= 4;
      a[0] = 0;
      a[1] = 1;
      b[0] = 0;
      b[1] = 1;
      c[0] = -1;
      c[1] = -1;
      a1d  = 1;
      b1d  = -3;
      c1d  = 2;
      aa[0] = 0;
      aa[1] = 4;    //changed from zero to 4
      bb[0][0] = 0;
      bb[0][1] = -4;
      bb[1][0] = -4;
      bb[1][1] = -8;
      out[10][0][0] = (aa[0] + bb[0][0]*in[0] + bb[1][0]*in[1]) * (a1d + in[2]*b1d + in[2]*in[2]*c1d);
      out[10][0][1] = (aa[1] + bb[0][1]*in[0] + bb[1][1]*in[1]) * (a1d + in[2]*b1d + in[2]*in[2]*c1d);
      out[10][0][2] = coeff * (a[0] + b[0]*in[0] + b[1]*in[1]) * (a[1] + c[0]*in[0] + c[1]*in[1]) * (b1d + 2*c1d*in[2]);



      coeff= 4;
      a[0] = 0;
      a[1] = 0;
      b[0] = 1;
      b[1] = 0;
      c[0] = 0;
      c[1] = 1;
      a1d  = 1;
      b1d  = -3;
      c1d  = 2;
      aa[0] = 0;
      aa[1] = 0;
      bb[0][0] = 0;
      bb[0][1] = 4;
      bb[1][0] = 4;
      bb[1][1] = 0;
      out[11][0][0] = (aa[0] + bb[0][0]*in[0] + bb[1][0]*in[1]) * (a1d + in[2]*b1d + in[2]*in[2]*c1d);
      out[11][0][1] = (aa[1] + bb[0][1]*in[0] + bb[1][1]*in[1]) * (a1d + in[2]*b1d + in[2]*in[2]*c1d);
      out[11][0][2] = coeff * (a[0] + b[0]*in[0] + b[1]*in[1]) * (a[1] + c[0]*in[0] + c[1]*in[1]) * (b1d + 2*c1d*in[2]);



      // upper triangle edges
      coeff= 4;
      a[0] = 0;
      a[1] = 1;
      b[0] = 1;
      b[1] = 0;
      c[0] = -1;
      c[1] = -1;
      a1d  = 0;
      b1d  = -1;
      c1d  = 2;
      aa[0] = 4;
      aa[1] = 0;
      bb[0][0] = -8;
      bb[0][1] = -4;
      bb[1][0] = -4;
      bb[1][1] = 0;
      out[12][0][0] = (aa[0] + bb[0][0]*in[0] + bb[1][0]*in[1]) * (a1d + in[2]*b1d + in[2]*in[2]*c1d);
      out[12][0][1] = (aa[1] + bb[0][1]*in[0] + bb[1][1]*in[1]) * (a1d + in[2]*b1d + in[2]*in[2]*c1d);
      out[12][0][2] = coeff * (a[0] + b[0]*in[0] + b[1]*in[1]) * (a[1] + c[0]*in[0] + c[1]*in[1]) * (b1d + 2*c1d*in[2]);



      coeff= 4;
      a[0] = 0;
      a[1] = 1;
      b[0] = 0;
      b[1] = 1;
      c[0] = -1;
      c[1] = -1;
      a1d  = 0;
      b1d  = -1;
      c1d  = 2;
      aa[0] = 0;
      aa[1] = 4;
      bb[0][0] = 0;
      bb[0][1] = -4;
      bb[1][0] = -4;
      bb[1][1] = -8;
      out[13][0][0] = (aa[0] + bb[0][0]*in[0] + bb[1][0]*in[1]) * (a1d + in[2]*b1d + in[2]*in[2]*c1d);
      out[13][0][1] = (aa[1] + bb[0][1]*in[0] + bb[1][1]*in[1]) * (a1d + in[2]*b1d + in[2]*in[2]*c1d);
      out[13][0][2] = coeff * (a[0] + b[0]*in[0] + b[1]*in[1]) * (a[1] + c[0]*in[0] + c[1]*in[1]) * (b1d + 2*c1d*in[2]);



      coeff= 4;
      a[0] = 0;
      a[1] = 0;
      b[0] = 1;
      b[1] = 0;
      c[0] = 0;
      c[1] = 1;
      a1d  = 0;
      b1d  = -1;
      c1d  = 2;
      aa[0] = 0;
      aa[1] = 0;
      bb[0][0] = 0;
      bb[0][1] = 4;
      bb[1][0] = 4;
      bb[1][1] = 0;
      out[14][0][0] = (aa[0] + bb[0][0]*in[0] + bb[1][0]*in[1]) * (a1d + in[2]*b1d + in[2]*in[2]*c1d);
      out[14][0][1] = (aa[1] + bb[0][1]*in[0] + bb[1][1]*in[1]) * (a1d + in[2]*b1d + in[2]*in[2]*c1d);
      out[14][0][2] = coeff * (a[0] + b[0]*in[0] + b[1]*in[1]) * (a[1] + c[0]*in[0] + c[1]*in[1]) * (b1d + 2*c1d*in[2]);



      // quadrilateral sides
      coeff= 4;
      a[0] = 0;
      a[1] = 1;
      b[0] = 1;
      b[1] = 0;
      c[0] = -1;
      c[1] = -1;
      a1d  = 0;
      b1d  = 4;
      c1d  = -4;
      aa[0] = 4;
      aa[1] = 0;
      bb[0][0] = -8;
      bb[0][1] = -4;
      bb[1][0] = -4;
      bb[1][1] = 0;
      out[15][0][0] = (aa[0] + bb[0][0]*in[0] + bb[1][0]*in[1]) * (a1d + in[2]*b1d + in[2]*in[2]*c1d);
      out[15][0][1] = (aa[1] + bb[0][1]*in[0] + bb[1][1]*in[1]) * (a1d + in[2]*b1d + in[2]*in[2]*c1d);
      out[15][0][2] = coeff * (a[0] + b[0]*in[0] + b[1]*in[1]) * (a[1] + c[0]*in[0] + c[1]*in[1]) * (b1d + 2*c1d*in[2]);





      coeff= 4;
      a[0] = 0;
      a[1] = 1;
      b[0] = 0;
      b[1] = 1;
      c[0] = -1;
      c[1] = -1;
      a1d  = 0;
      b1d  = 4;
      c1d  = -4;
      aa[0] = 0;
      aa[1] = 4;
      bb[0][0] = 0;
      bb[0][1] = -4;
      bb[1][0] = -4;
      bb[1][1] = -8;
      out[16][0][0] = (aa[0] + bb[0][0]*in[0] + bb[1][0]*in[1]) * (a1d + in[2]*b1d + in[2]*in[2]*c1d);
      out[16][0][1] = (aa[1] + bb[0][1]*in[0] + bb[1][1]*in[1]) * (a1d + in[2]*b1d + in[2]*in[2]*c1d);
      out[16][0][2] = coeff * (a[0] + b[0]*in[0] + b[1]*in[1]) * (a[1] + c[0]*in[0] + c[1]*in[1]) * (b1d + 2*c1d*in[2]);




      coeff= 4;
      a[0] = 0;
      a[1] = 0;
      b[0] = 1;
      b[1] = 0;
      c[0] = 0;
      c[1] = 1;
      a1d  = 0;
      b1d  = 4;
      c1d  = -4;
      aa[0] = 0;
      aa[1] = 0;
      bb[0][0] = 0;
      bb[0][1] = 4;
      bb[1][0] = 4;
      bb[1][1] = 0;
      out[17][0][0] = (aa[0] + bb[0][0]*in[0] + bb[1][0]*in[1]) * (a1d + in[2]*b1d + in[2]*in[2]*c1d);
      out[17][0][1] = (aa[1] + bb[0][1]*in[0] + bb[1][1]*in[1]) * (a1d + in[2]*b1d + in[2]*in[2]*c1d);
      out[17][0][2] = coeff * (a[0] + b[0]*in[0] + b[1]*in[1]) * (a[1] + c[0]*in[0] + c[1]*in[1]) * (b1d + 2*c1d*in[2]);

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
        out.resize(size());

        auto const direction = std::distance(order.begin(), std::find(order.begin(), order.end(), 1));
        switch (direction) {
        case 0:
          out[0] = (-3 + 4*(in[0] + in[1])) * (1 - 3*in[2] + 2*in[2]*in[2]);
          out[1] = (-1 + 4*in[0]) * (1 - 3*in[2] + 2*in[2]*in[2]);
          out[2] = 0;
          out[3] = (-3 + 4*(in[0] + in[1])) * (-in[2] + 2*in[2]*in[2]);
          out[4] = (-1 + 4*in[0]) * (-in[2] + 2*in[2]*in[2]);
          out[5] = 0;
          out[6] = (-3 + 4*(in[0] + in[1])) * 4*(in[2] - in[2]*in[2]);
          out[7] = (-1 + 4*in[0]) * 4*(in[2] - in[2]*in[2]);
          out[8] = 0;
          out[9] = (4 - 8*in[0] - 4*in[1]) * (1 - 3*in[2] + 2*in[2]*in[2]);
          out[10] = -4*in[1] * (1 - 3*in[2] + 2*in[2]*in[2]);
          out[11] = 4*in[1] * (1 - 3*in[2] + 2*in[2]*in[2]);
          out[12] = (4 - 8*in[0] - 4*in[1]) * (-in[2] + 2*in[2]*in[2]);
          out[13] = -4*in[1] * (-in[2] + 2*in[2]*in[2]);
          out[14] = 4*in[1] * (-in[2] + 2*in[2]*in[2]);
          out[15] = (4 - 8*in[0] - 4*in[1]) * (4*in[2] - 4*in[2]*in[2]);
          out[16] = -4*in[1] * (4*in[2] - 4*in[2]*in[2]);
          out[17] = 4*in[1] * (4*in[2] - 4*in[2]*in[2]);
          break;

        case 1:
          out[0] = (-3 + 4*(in[0] + in[1])) * (1 - 3*in[2] + 2*in[2]*in[2]);
          out[1] = 0;
          out[2] = (-1 + 4*in[1]) * (1 - 3*in[2] + 2*in[2]*in[2]);
          out[3] = (-3 + 4*(in[0] + in[1])) * (-in[2] + 2*in[2]*in[2]);
          out[4] = 0;
          out[5] = (-1 + 4*in[1]) * (-in[2] + 2*in[2]*in[2]);
          out[6] = (-3 + 4*(in[0] + in[1])) * 4*(in[2] - in[2]*in[2]);
          out[7] = 0;
          out[8] = (-1 + 4*in[1]) * 4*(in[2] - in[2]*in[2]);
          out[9] = -4*in[0] * (1 - 3*in[2] + 2*in[2]*in[2]);
          out[10] = (4 - 4*in[0] - 8*in[1]) * (1 - 3*in[2] + 2*in[2]*in[2]);
          out[11] = 4*in[0] * (1 - 3*in[2] + 2*in[2]*in[2]);
          out[12] = -4*in[0] * (-in[2] + 2*in[2]*in[2]);
          out[13] = (4 - 4*in[0] - 8*in[1]) * (-in[2] + 2*in[2]*in[2]);
          out[14] = 4*in[0] * (-in[2] + 2*in[2]*in[2]);
          out[15] = -4*in[0] * 4*(in[2] - in[2]*in[2]);
          out[16] = (4 - 4*in[0] - 8*in[1]) * (4*in[2] - 4*in[2]*in[2]);
          out[17] = 4*in[0] * (4*in[2] - 4*in[2]*in[2]);
          break;

        case 2:
          out[0] = 2 * (1 - in[0] - in[1]) * (0.5 - in[0] - in[1]) * (-3 + 4*in[2]);
          out[1] = 2 * in[0] * (-0.5 + in[0]) * (-3 + 4*in[2]);
          out[2] = 2 * in[1] * (-0.5 + in[1]) * (-3 + 4*in[2]);
          out[3] = 2 * (1 - in[0] - in[1]) * (0.5 - in[0] - in[1]) * (-1 + 4*in[2]);
          out[4] = 2 * in[0] * (-0.5 + in[0]) * (-1 + 4*in[2]);
          out[5] = 2 * in[1] * (-0.5 + in[1]) * (-1 + 4*in[2]);
          out[6] = 2 * (1 - in[0] - in[1]) * (0.5 - in[0] - in[1]) * (4 - 8*in[2]);
          out[7] = 2 * in[0] * (-0.5 + in[0]) * (4 - 8*in[2]);
          out[8] = 2*in[1] * (-0.5 + in[1]) * (4 - 8*in[2]);
          out[9] =  4*in[0] * (1 - in[0] - in[1]) * (-3 + 4*in[2]);
          out[10] = 4*in[1] * (1 - in[0] - in[1]) * (-3 + 4*in[2]);
          out[11] = 4*in[0]*in[1] * (-3 + 4*in[2]);
          out[12] = 4 * in[0] * (1 - in[0] - in[1]) * (-1 + 4*in[2]);
          out[13] = 4 * in[1] * (1 - in[0] - in[1]) * (-1 + 4*in[2]);
          out[14] = 4*in[0]*in[1] * (-1 + 4*in[2]);
          out[15] = 4 * in[0] * (1 - in[0] - in[1]) * (4 - 8*in[2]);
          out[16] = 4 * in[1] * (1 - in[0] - in[1]) * (4 - 8*in[2]);
          out[17] = 4*in[0]*in[1] * (4 - 8*in[2]);
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
