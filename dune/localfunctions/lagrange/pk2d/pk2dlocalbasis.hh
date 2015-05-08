// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PK2DLOCALBASIS_HH
#define DUNE_PK2DLOCALBASIS_HH

#include <dune/common/fmatrix.hh>

#include <dune/localfunctions/common/localbasis.hh>

namespace Dune
{
  /**@ingroup LocalBasisImplementation
         \brief Lagrange shape functions of arbitrary order on the reference triangle.

         Lagrange shape functions of arbitrary order have the property that
         \f$\hat\phi^i(x_j) = \delta_{i,j}\f$ for certain points \f$x_j\f$.

         \tparam D Type to represent the field in the domain.
         \tparam R Type to represent the field in the range.
         \tparam k Polynomial order.

         \nosubgrouping
   */
  template<class D, class R, unsigned int k>
  class Pk2DLocalBasis
  {
  public:

    /** \brief Export the number of degrees of freedom */
    enum {N = (k+1)*(k+2)/2};

    /** \brief Export the element order
       OS: Surprising that we need to export this both statically and dynamically!
     */
    enum {O = k};

    typedef LocalBasisTraits<D,2,Dune::FieldVector<D,2>,R,1,Dune::FieldVector<R,1>,
        Dune::FieldMatrix<R,1,2>, 2 > Traits;

    //! \brief Standard constructor
    Pk2DLocalBasis ()
    {
      for (unsigned int i=0; i<=k; i++)
        pos[i] = (1.0*i)/std::max(k,(unsigned int)1);
    }

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
      // specialization for k==0, not clear whether that is needed
      if (k==0) {
        out[0] = 1;
        return;
      }

      int n=0;
      for (unsigned int j=0; j<=k; j++)
        for (unsigned int i=0; i<=k-j; i++)
        {
          out[n] = 1.0;
          for (unsigned int alpha=0; alpha<i; alpha++)
            out[n] *= (x[0]-pos[alpha])/(pos[i]-pos[alpha]);
          for (unsigned int beta=0; beta<j; beta++)
            out[n] *= (x[1]-pos[beta])/(pos[j]-pos[beta]);
          for (unsigned int gamma=i+j+1; gamma<=k; gamma++)
            out[n] *= (pos[gamma]-x[0]-x[1])/(pos[gamma]-pos[i]-pos[j]);
          n++;
        }
    }

    //! \brief Evaluate Jacobian of all shape functions
    inline void
    evaluateJacobian (const typename Traits::DomainType& x,       // position
                      std::vector<typename Traits::JacobianType>& out) const                        // return value
    {
      out.resize(N);

      // specialization for k==0, not clear whether that is needed
      if (k==0) {
        out[0][0][0] = 0; out[0][0][1] = 0;
        return;
      }

      int n=0;
      for (unsigned int j=0; j<=k; j++)
        for (unsigned int i=0; i<=k-j; i++)
        {
          // x_0 derivative
          out[n][0][0] = 0.0;
          R factor=1.0;
          for (unsigned int beta=0; beta<j; beta++)
            factor *= (x[1]-pos[beta])/(pos[j]-pos[beta]);
          for (unsigned int a=0; a<i; a++)
          {
            R product=factor;
            for (unsigned int alpha=0; alpha<i; alpha++)
              if (alpha==a)
                product *= 1.0/(pos[i]-pos[alpha]);
              else
                product *= (x[0]-pos[alpha])/(pos[i]-pos[alpha]);
            for (unsigned int gamma=i+j+1; gamma<=k; gamma++)
              product *= (pos[gamma]-x[0]-x[1])/(pos[gamma]-pos[i]-pos[j]);
            out[n][0][0] += product;
          }
          for (unsigned int c=i+j+1; c<=k; c++)
          {
            R product=factor;
            for (unsigned int alpha=0; alpha<i; alpha++)
              product *= (x[0]-pos[alpha])/(pos[i]-pos[alpha]);
            for (unsigned int gamma=i+j+1; gamma<=k; gamma++)
              if (gamma==c)
                product *= -1.0/(pos[gamma]-pos[i]-pos[j]);
              else
                product *= (pos[gamma]-x[0]-x[1])/(pos[gamma]-pos[i]-pos[j]);
            out[n][0][0] += product;
          }

          // x_1 derivative
          out[n][0][1] = 0.0;
          factor = 1.0;
          for (unsigned int alpha=0; alpha<i; alpha++)
            factor *= (x[0]-pos[alpha])/(pos[i]-pos[alpha]);
          for (unsigned int b=0; b<j; b++)
          {
            R product=factor;
            for (unsigned int beta=0; beta<j; beta++)
              if (beta==b)
                product *= 1.0/(pos[j]-pos[beta]);
              else
                product *= (x[1]-pos[beta])/(pos[j]-pos[beta]);
            for (unsigned int gamma=i+j+1; gamma<=k; gamma++)
              product *= (pos[gamma]-x[0]-x[1])/(pos[gamma]-pos[i]-pos[j]);
            out[n][0][1] += product;
          }
          for (unsigned int c=i+j+1; c<=k; c++)
          {
            R product=factor;
            for (unsigned int beta=0; beta<j; beta++)
              product *= (x[1]-pos[beta])/(pos[j]-pos[beta]);
            for (unsigned int gamma=i+j+1; gamma<=k; gamma++)
              if (gamma==c)
                product *= -1.0/(pos[gamma]-pos[i]-pos[j]);
              else
                product *= (pos[gamma]-x[0]-x[1])/(pos[gamma]-pos[i]-pos[j]);
            out[n][0][1] += product;
          }

          n++;
        }

    }

    //! \brief Evaluate higher derivatives of all shape functions
    template<unsigned int dOrder> //order of derivative
    inline void evaluate(const std::array<int,dOrder>& directions, //direction of derivative
                         const typename Traits::DomainType& in,  //position
                         std::vector<typename Traits::RangeType>& out) const //return value
    {
      out.resize(N);

      if (dOrder > Traits::diffOrder)
        DUNE_THROW(NotImplemented, "Desired derivative order is not implemented");

      if (dOrder==0)
        evaluateFunction(in, out);
      else if (dOrder==1)
      {
        int n=0;
        for (unsigned int j=0; j<=k; j++)
          for (unsigned int i=0; i<=k-j; i++, n++)
          {
            out[n] = 0.0;
            for (unsigned int no1=0; no1 < k; no1++)
            {
              R factor = lagrangianFactorDerivative(directions[0], no1, i, j, in);
              for (unsigned int no2=0; no2 < k; no2++)
                if (no1 != no2)
                  factor *= lagrangianFactor(no2, i, j, in);

              out[n] += factor;
            }
          }
      }
      else if (dOrder==2)
      {
        // specialization for k<2, not clear whether that is needed
        if (k<2) {
          std::fill(out.begin(), out.end(), 0.0);
        return;
      }

      //f = prod_{i} f_i -> dxa dxb f = sum_{i} {dxa f_i sum_{k \neq i} dxb f_k prod_{l \neq k,i} f_l
      int n=0;
      for (unsigned int j=0; j<=k; j++)
        for (unsigned int i=0; i<=k-j; i++, n++)
        {
          R res = 0.0;

          for (unsigned int no1=0; no1 < k; no1++)
          {
            R factor1 = lagrangianFactorDerivative(directions[0], no1, i, j, in);
            for (unsigned int no2=0; no2 < k; no2++)
            {
              if (no1 == no2)
                continue;
              R factor2 = factor1*lagrangianFactorDerivative(directions[1], no2, i, j, in);
              for (unsigned int no3=0; no3 < k; no3++)
              {
                if (no3 == no1 || no3 == no2)
                  continue;
                factor2 *= lagrangianFactor(no3, i, j, in);
              }
              res += factor2;
            }

          }
          out[n] = res;
        }
      }
    }

    //! \brief Polynomial order of the shape functions
    unsigned int order () const
    {
      return k;
    }

  private:
  /** \brief Returns a single Lagrangian factor of l_ij evaluated at x */
  typename Traits::RangeType lagrangianFactor(const int no, const int i, const int j, const typename Traits::DomainType& x) const
  {
    if ( no < i)
      return (x[0]-pos[no])/(pos[i]-pos[no]);
    if (no < i+j)
      return (x[1]-pos[no-i])/(pos[j]-pos[no-i]);
    return (pos[no+1]-x[0]-x[1])/(pos[no+1]-pos[i]-pos[j]);
  }

  /** \brief Returns the derivative of a single Lagrangian factor of l_ij evaluated at x
   * \param direction Derive in x-direction if this is 0, otherwise derive in y direction
   */
  typename Traits::RangeType lagrangianFactorDerivative(const int direction, const int no, const int i, const int j, const typename Traits::DomainType& x) const
  {
    if ( no < i)
      return (direction == 0) ? 1.0/(pos[i]-pos[no]) : 0;

    if (no < i+j)
      return (direction == 0) ? 0: 1.0/(pos[j]-pos[no-i]);

    return -1.0/(pos[no+1]-pos[i]-pos[j]);
  }

    R pos[k+1]; // positions on the interval
  };

}
#endif
