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
        Dune::FieldMatrix<R,1,2> > Traits;

    //! \brief Standard constructor
    Pk2DLocalBasis ()
    {
      for (unsigned int i=0; i<=k; i++)
        pos[i] = (1.0*i)/k;
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

    //! \brief Polynomial order of the shape functions
    unsigned int order () const
    {
      return k;
    }

  private:
    R pos[k+1]; // positions on the interval
  };

}
#endif
