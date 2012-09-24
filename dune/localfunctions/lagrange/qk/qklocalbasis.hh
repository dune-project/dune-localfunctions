// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_LOCALFUNCTIONS_QKLOCALBASIS_HH
#define DUNE_LOCALFUNCTIONS_QKLOCALBASIS_HH

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/misc.hh>

#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/localbasis.hh>
#include <dune/localfunctions/common/localfiniteelementtraits.hh>


namespace Dune
{
  /**@ingroup LocalBasisImplementation
     \brief Lagrange shape functions of order k on the reference cube.

     Also known as \f$Q^k\f$.

     \tparam D Type to represent the field in the domain.
     \tparam R Type to represent the field in the range.
     \tparam k Polynomial degree
     \tparam d Dimension of the cube

     \nosubgrouping
   */
  template<class D, class R, int k, int d>
  class QkLocalBasis
  {
    enum { n = Power_m_p<k+1,d>::power };

    // ith Lagrange polynomial of degree k in one dimension
    static R p (int i, D x)
    {
      R result(1.0);
      for (int j=0; j<=k; j++)
        if (j!=i) result *= (k*x-j)/(i-j);
      return result;
    }

    // derivative of ith Lagrange polynomial of degree k in one dimension
    static R dp (int i, D x)
    {
      R result(0.0);

      for (int j=0; j<=k; j++)
        if (j!=i)
        {
          R prod( (k*1.0)/(i-j) );
          for (int l=0; l<=k; l++)
            if (l!=i && l!=j)
              prod *= (k*x-l)/(i-l);
          result += prod;
        }
      return result;
    }

    // Return i as a d-digit number in the (k+1)-nary system
    static Dune::FieldVector<int,d> multiindex (int i)
    {
      Dune::FieldVector<int,d> alpha;
      for (int j=0; j<d; j++)
      {
        alpha[j] = i % (k+1);
        i = i/(k+1);
      }
      return alpha;
    }

  public:
    typedef LocalBasisTraits<D,d,Dune::FieldVector<D,d>,R,1,Dune::FieldVector<R,1>,Dune::FieldMatrix<R,1,d> > Traits;

    //! \brief number of shape functions
    unsigned int size () const
    {
      return Power_m_p<k+1,d>::power;
    }

    //! \brief Evaluate all shape functions
    inline void evaluateFunction (const typename Traits::DomainType& in,
                                  std::vector<typename Traits::RangeType>& out) const
    {
      out.resize(size());
      for (size_t i=0; i<size(); i++)
      {
        // convert index i to multiindex
        Dune::FieldVector<int,d> alpha(multiindex(i));

        // initialize product
        out[i] = 1.0;

        // dimension by dimension
        for (int j=0; j<d; j++)
          out[i] *= p(alpha[j],in[j]);
      }
    }

    /** \brief Evaluate Jacobian of all shape functions
     * \param in position where to evaluate
     * \param out The return value
     */
    inline void
    evaluateJacobian (const typename Traits::DomainType& in,
                      std::vector<typename Traits::JacobianType>& out) const
    {
      out.resize(size());

      // Loop over all shape functions
      for (size_t i=0; i<size(); i++)
      {
        // convert index i to multiindex
        Dune::FieldVector<int,d> alpha(multiindex(i));

        // Loop over all coordinate directions
        for (int j=0; j<d; j++)
        {
          // Initialize: the overall expression is a product
          // if j-th bit of i is set to -1, else 1
          out[i][0][j] = dp(alpha[j],in[j]);

          // rest of the product
          for (int l=0; l<d; l++)
            if (l!=j)
              out[i][0][j] *= p(alpha[l],in[l]);
        }
      }
    }

    //! \brief Polynomial order of the shape functions
    unsigned int order () const
    {
      return k;
    }
  };
}

#endif
