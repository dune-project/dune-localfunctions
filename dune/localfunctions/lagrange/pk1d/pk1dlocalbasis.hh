// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_Pk1DLOCALBASIS_HH
#define DUNE_Pk1DLOCALBASIS_HH

#include <numeric>

#include <dune/common/fmatrix.hh>

#include <dune/localfunctions/common/localbasis.hh>

namespace Dune
{
  /**@ingroup LocalBasisImplementation
     \brief Lagrange shape functions of arbitrary order on the 1D reference triangle.

     Lagrange shape functions of arbitrary order have the property that
     \f$\hat\phi^i(x_j) = \delta_{i,j}\f$ for certain points \f$x_j\f$.

     \tparam D Type to represent the field in the domain.
     \tparam R Type to represent the field in the range.
     \tparam k Polynomial order.

     \nosubgrouping
   */
  template<class D, class R, unsigned int k>
  class Pk1DLocalBasis
  {
  public:

    /** \brief Export the number of degrees of freedom */
    enum {N = k+1};

    /** \brief Export the element order */
    enum {O = k};

    typedef LocalBasisTraits<D,
        1,
        Dune::FieldVector<D,1>,
        R,
        1,
        Dune::FieldVector<R,1>,
        Dune::FieldMatrix<R,1,1>,
        1> Traits;

    //! \brief Standard constructor
    Pk1DLocalBasis ()
    {
      for (unsigned int i=0; i<=k; i++)
        pos[i] = (1.0*i)/std::max(k,(unsigned int)1);
    }

    //! \brief number of shape functions
    constexpr unsigned int size () const
    {
      return N;
    }

    //! \brief Evaluate all shape functions
    inline void evaluateFunction (const typename Traits::DomainType& x,
                                  std::vector<typename Traits::RangeType>& out) const
    {
      out.resize(N);

      for (unsigned int i=0; i<N; i++)
      {
        out[i] = 1.0;
        for (unsigned int alpha=0; alpha<i; alpha++)
          out[i] *= (x[0]-pos[alpha])/(pos[i]-pos[alpha]);
        for (unsigned int gamma=i+1; gamma<=k; gamma++)
          out[i] *= (x[0]-pos[gamma])/(pos[i]-pos[gamma]);
      }
    }


    //! \brief Evaluate Jacobian of all shape functions.
    inline void
    evaluateJacobian (const typename Traits::DomainType& x,                                      // position
                      std::vector<typename Traits::JacobianType>& out) const                             // return value
    {
      out.resize(N);
      evaluateJacobianTemplate(x, [&out](unsigned int i)->R& { return out[i][0][0]; });
    }


    /** \brief Evaluate partial derivatives of any order of all shape functions
     * \param order Order of the partial derivatives, in the classic multi-index notation
     * \param in Position where to evaluate the derivatives
     * \param[out] out Return value: the desired partial derivatives
     */
    void partial(const std::array<unsigned int,1>& order,
                 const typename Traits::DomainType& in,
                 std::vector<typename Traits::RangeType>& out) const
    {
      auto totalOrder = order[0];
      if (totalOrder == 0) {
        evaluateFunction(in, out);
      } else if (totalOrder == 1) {
        out.resize(N);
        evaluateJacobianTemplate(in, [&out](unsigned int i)->R& { return out[i]; });
      } else {
        DUNE_THROW(NotImplemented, "Desired derivative order is not implemented");
      }
    }

    //! \brief Evaluate higher derivatives of all shape functions
    template<std::size_t dOrder> //order of derivative
    inline void evaluate(const std::array<int,dOrder>& /*directions*/,                              // direction of derivative
                         const typename Traits::DomainType& in,                                 // position
                         std::vector<typename Traits::RangeType>& out) const    // return value
    {
      std::array<unsigned int,1> order;
      order[0] = dOrder;
      partial(order, in, out);
    }

    //! \brief Polynomial order of the shape functions
    constexpr unsigned int order () const
    {
      return k;
    }

  private:

    // Evaluate Jacobian of all shape functions. Use an assigner argument as
    // return type, i.e. out(i) must be mutable reference to i'th jacobian entry.
    template<class Assigner>
    inline void evaluateJacobianTemplate (const typename Traits::DomainType& x, // position
                                          Assigner out) const                                                                           // return value
    {
      // expects the out-vector to be resized to size == N
      for (unsigned int i=0; i<=k; i++) {

        // x_0 derivative
        out(i) = 0.0;
        R factor=1.0;
        for (unsigned int a=0; a<i; a++)
        {
          R product=factor;
          for (unsigned int alpha=0; alpha<i; alpha++)
            product *= (alpha==a) ? 1.0/(pos[i]-pos[alpha])
                       : (x[0]-pos[alpha])/(pos[i]-pos[alpha]);
          for (unsigned int gamma=i+1; gamma<=k; gamma++)
            product *= (pos[gamma]-x[0])/(pos[gamma]-pos[i]);
          out(i) += product;
        }
        for (unsigned int c=i+1; c<=k; c++)
        {
          R product=factor;
          for (unsigned int alpha=0; alpha<i; alpha++)
            product *= (x[0]-pos[alpha])/(pos[i]-pos[alpha]);
          for (unsigned int gamma=i+1; gamma<=k; gamma++)
            product *= (gamma==c) ? -1.0/(pos[gamma]-pos[i])
                       : (pos[gamma]-x[0])/(pos[gamma]-pos[i]);
          out(i) += product;
        }
      }
    }

  private:
    R pos[k+1];     // positions on the interval
  };

}
#endif
