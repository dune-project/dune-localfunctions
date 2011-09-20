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


      // Compute function values: they are products of 1d Lagrange function values
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

      // Evaluate the 1d Lagrange functions and their derivatives
      array<array<R,3>, dim> X, DX;

      for (size_t i=0; i<dim; i++) {
        X[i][0] =  R(2)*in[i]*in[i] - R(3)*in[i]+R(1);
        X[i][1] = -R(4)*in[i]*in[i] + R(4)*in[i];
        X[i][2] =  R(2)*in[i]*in[i] -   in[i];

        DX[i][0] =  R(4)*in[i] - R(3);
        DX[i][1] = -R(8)*in[i] + R(4);
        DX[i][2] =  R(4)*in[i] - R(1);
      }


      // Compute the derivatives by deriving the products of 1d Lagrange functions
      for (size_t i=0; i<out.size(); i++) {

        // Computing the j-th partial derivative
        for (int j=0; j<dim; j++) {

          out[i][0][j] = 1;

          // Loop over the 'dim' terms in the product rule
          size_t ternary = i;
          for (int k=0; k<dim; k++) {

            int digit = ternary%3;
            ternary /= 3;

            out[i][0][j] *= (k==j) ? DX[k][digit] : X[k][digit];

          }

        }

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
