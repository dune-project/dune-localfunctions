// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PK3DLOCALBASIS_HH
#define DUNE_PK3DLOCALBASIS_HH

#include <numeric>

#include <dune/common/fmatrix.hh>

#include <dune/localfunctions/common/localbasis.hh>

namespace Dune
{
  /**@ingroup LocalBasisImplementation
     \brief Lagrange shape functions of arbitrary order on the reference tetrahedron.

     Lagrange shape functions of arbitrary order have the property that
     \f$\hat\phi^i(x_j) = \delta_{i,j}\f$ for certain points \f$x_j\f$.

     \tparam D Type to represent the field in the domain.
     \tparam R Type to represent the field in the range.
     \tparam k Polynomial order.

     \nosubgrouping
   */
  template<class D, class R, unsigned int k>
  class Pk3DLocalBasis
  {
  public:
    enum {N = (k+1)*(k+2)*(k+3)/6};
    enum {O = k};

    typedef LocalBasisTraits<D,3,Dune::FieldVector<D,3>,R,1,Dune::FieldVector<R,1>,
        Dune::FieldMatrix<R,1,3> > Traits;

    //! \brief Standard constructor
    Pk3DLocalBasis () {}

    //! \brief number of shape functions
    unsigned int size () const
    {
      return N;
    }

    //! \brief Evaluate all shape functions
    inline void evaluateFunction (const typename Traits::DomainType& x,
                                  std::vector<typename Traits::RangeType>& out) const
    {
      out.resize(N);
      typename Traits::DomainType kx = x;
      kx *= k;
      unsigned int n = 0;
      unsigned int i[4];
      R factor[4];
      for (i[2] = 0; i[2] <= k; ++i[2])
      {
        factor[2] = 1.0;
        for (unsigned int j = 0; j < i[2]; ++j)
          factor[2] *= (kx[2]-j) / (i[2]-j);
        for (i[1] = 0; i[1] <= k - i[2]; ++i[1])
        {
          factor[1] = 1.0;
          for (unsigned int j = 0; j < i[1]; ++j)
            factor[1] *= (kx[1]-j) / (i[1]-j);
          for (i[0] = 0; i[0] <= k - i[1] - i[2]; ++i[0])
          {
            factor[0] = 1.0;
            for (unsigned int j = 0; j < i[0]; ++j)
              factor[0] *= (kx[0]-j) / (i[0]-j);
            i[3] = k - i[0] - i[1] - i[2];
            D kx3 = k - kx[0] - kx[1] - kx[2];
            factor[3] = 1.0;
            for (unsigned int j = 0; j < i[3]; ++j)
              factor[3] *= (kx3-j) / (i[3]-j);
            out[n++] = factor[0] * factor[1] * factor[2] * factor[3];
          }
        }
      }
    }

    //! \brief Evaluate Jacobian of all shape functions
    inline void
    evaluateJacobian (const typename Traits::DomainType& x,         // position
                      std::vector<typename Traits::JacobianType>& out) const      // return value
    {
      out.resize(N);
      typename Traits::DomainType kx = x;
      kx *= k;
      unsigned int n = 0;
      unsigned int i[4];
      R factor[4];
      for (i[2] = 0; i[2] <= k; ++i[2])
      {
        factor[2] = 1.0;
        for (unsigned int j = 0; j < i[2]; ++j)
          factor[2] *= (kx[2]-j) / (i[2]-j);
        for (i[1] = 0; i[1] <= k - i[2]; ++i[1])
        {
          factor[1] = 1.0;
          for (unsigned int j = 0; j < i[1]; ++j)
            factor[1] *= (kx[1]-j) / (i[1]-j);
          for (i[0] = 0; i[0] <= k - i[1] - i[2]; ++i[0])
          {
            factor[0] = 1.0;
            for (unsigned int j = 0; j < i[0]; ++j)
              factor[0] *= (kx[0]-j) / (i[0]-j);
            i[3] = k - i[0] - i[1] - i[2];
            D kx3 = k - kx[0] - kx[1] - kx[2];
            R sum3 = 0.0;
            factor[3] = 1.0;
            for (unsigned int j = 0; j < i[3]; ++j)
              factor[3] /= i[3] - j;
            R prod_all = factor[0] * factor[1] * factor[2] * factor[3];
            for (unsigned int j = 0; j < i[3]; ++j)
            {
              R prod = prod_all;
              for (unsigned int l = 0; l < i[3]; ++l)
                if (j == l)
                  prod *= -R(k);
                else
                  prod *= kx3 - l;
              sum3 += prod;
            }
            for (unsigned int j = 0; j < i[3]; ++j)
              factor[3] *= kx3 - j;
            for (unsigned int m = 0; m < 3; ++m)
            {
              out[n][0][m] = sum3;
              for (unsigned int j = 0; j < i[m]; ++j)
              {
                R prod = factor[3];
                for (unsigned int p = 0; p < 3; ++p)
                {
                  if (m == p)
                    for (unsigned int l = 0; l < i[p]; ++l)
                      if (j == l)
                        prod *= R(k) / (i[p]-l);
                      else
                        prod *= (kx[p]-l) / (i[p]-l);
                  else
                    prod *= factor[p];
                }
                out[n][0][m] += prod;
              }
            }
            n++;
          }
        }
      }
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
      } else {
        DUNE_THROW(NotImplemented, "Desired derivative order is not implemented");
      }
    }

    //! \brief Polynomial order of the shape functions
    unsigned int order () const
    {
      return k;
    }
  };


  //Specialization for k=0
  template<class D, class R>
  class Pk3DLocalBasis<D,R,0>
  {
  public:
    typedef LocalBasisTraits<D,3,Dune::FieldVector<D,3>,R,1,Dune::FieldVector<R,1>,
        Dune::FieldMatrix<R,1,3> > Traits;

    /** \brief Export the number of degrees of freedom */
    enum {N = 1};

    /** \brief Export the element order */
    enum {O = 0};

    unsigned int size () const
    {
      return 1;
    }

    inline void evaluateFunction (const typename Traits::DomainType& in,
                                  std::vector<typename Traits::RangeType>& out) const
    {
      out.resize(1);
      out[0] = 1;
    }

    // evaluate derivative of a single component
    inline void
    evaluateJacobian (const typename Traits::DomainType& in,         // position
                      std::vector<typename Traits::JacobianType>& out) const      // return value
    {
      out.resize(1);
      out[0][0][0] = 0;
      out[0][0][1] = 0;
      out[0][0][2] = 0;
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
      } else {
        out.resize(N);
        out[0] = 0;
      }
    }

    // local interpolation of a function
    template<typename E, typename F, typename C>
    void interpolate (const E& e, const F& f, std::vector<C>& out) const
    {
      typename Traits::DomainType x;
      typename Traits::RangeType y;
      x[0] = 1.0/4.0;
      x[1] = 1.0/4.0;
      x[2] = 1.0/4.0;
      f.eval_local(e,x,y);
      out[0] = y;
    }

    unsigned int order () const
    {
      return 0;
    }
  };
}
#endif
