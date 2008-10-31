// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PK2DLOCALBASIS_HH
#define DUNE_PK2DLOCALBASIS_HH

#include "localbasis.hh"

/**@ingroup LocalBasisImplementation
   \brief Lagrange shape functions of arbitrary order on the reference triangle.

   Lagrange shape functions of arbitrary order have the property that
   \f$\hat\phi^i(x_j) = \delta_{i,j}\f$ for certain points \f$x_j\f$.

   - <tt>D</tt>: Type to represent the field in the domain.
   - <tt>R</tt>: Type to represent the field in the range.
   - <tt>k</tt>: Polynomial order.

   \nosubgrouping
 */
template<class D, class R, unsigned int k>
class Pk2DLocalBasis :
  public C1LocalBasisInterface<
      C1LocalBasisTraits<D,2,Dune::FieldVector<D,2>,R,1,Dune::FieldVector<R,1>,
          Dune::FieldVector<Dune::FieldVector<R,2>,1> >,
      Pk2DLocalBasis<D,R,k> >
{
  enum {N = (k+1)*(k+2)/2};
public:
  typedef C1LocalBasisTraits<D,2,Dune::FieldVector<D,2>,R,1,Dune::FieldVector<R,1>,
      Dune::FieldVector<Dune::FieldVector<R,2>,1> > Traits;

  //! \brief Standard constructor
  Pk2DLocalBasis ()
  {
    for (int i=0; i<=k; i++) pos[i] = (1.0*i)/k;
  }

  //! \brief number of shape functions
  unsigned int size () const
  //! \brief number of shape functions
  {
    return N;
  }

  //! \brief Evaluate all shape functions
  inline void evaluateFunction (const typename Traits::DomainType& x,
                                std::vector<typename Traits::RangeType>& out) const
  {
    out.resize(N);
    int n=0;
    for (int j=0; j<=k; j++)
      for (int i=0; i<=k-j; i++)
      {
        out[n] = 1.0;
        for (int alpha=0; alpha<i; alpha++)
          out[n] *= (x[0]-pos[alpha])/(pos[i]-pos[alpha]);
        for (int beta=0; beta<j; beta++)
          out[n] *= (x[1]-pos[beta])/(pos[j]-pos[beta]);
        for (int gamma=i+j+1; gamma<=k; gamma++)
          out[n] *= (pos[gamma]-x[0]-x[1])/(pos[gamma]-pos[i]-pos[j]);
        n++;
      }
  }

  //! \brief Evaluate Jacobian of all shape functions
  inline void
  evaluateJacobian (const typename Traits::DomainType& x,         // position
                    std::vector<typename Traits::JacobianType>& out) const // return value
  {
    out.resize(N);
    int n=0;
    for (int j=0; j<=k; j++)
      for (int i=0; i<=k-j; i++)
      {
        // x_0 derivative
        out[n][0][0] = 0.0;
        R factor=1.0;
        for (int beta=0; beta<j; beta++)
          factor *= (x[1]-pos[beta])/(pos[j]-pos[beta]);
        for (int a=0; a<i; a++)
        {
          R product=factor;
          for (int alpha=0; alpha<i; alpha++)
            if (alpha==a)
              product *= 1.0/(pos[i]-pos[alpha]);
            else
              product *= (x[0]-pos[alpha])/(pos[i]-pos[alpha]);
          for (int gamma=i+j+1; gamma<=k; gamma++)
            product *= (pos[gamma]-x[0]-x[1])/(pos[gamma]-pos[i]-pos[j]);
          out[n][0][0] += product;
        }
        for (int c=i+j+1; c<=k; c++)
        {
          R product=factor;
          for (int alpha=0; alpha<i; alpha++)
            product *= (x[0]-pos[alpha])/(pos[i]-pos[alpha]);
          for (int gamma=i+j+1; gamma<=k; gamma++)
            if (gamma==c)
              product *= -1.0/(pos[gamma]-pos[i]-pos[j]);
            else
              product *= (pos[gamma]-x[0]-x[1])/(pos[gamma]-pos[i]-pos[j]);
          out[n][0][0] += product;
        }

        // x_1 derivative
        out[n][0][1] = 0.0;
        factor = 1.0;
        for (int alpha=0; alpha<i; alpha++)
          factor *= (x[0]-pos[alpha])/(pos[i]-pos[alpha]);
        for (int b=0; b<j; b++)
        {
          R product=factor;
          for (int beta=0; beta<j; beta++)
            if (beta==b)
              product *= 1.0/(pos[j]-pos[beta]);
            else
              product *= (x[1]-pos[beta])/(pos[j]-pos[beta]);
          for (int gamma=i+j+1; gamma<=k; gamma++)
            product *= (pos[gamma]-x[0]-x[1])/(pos[gamma]-pos[i]-pos[j]);
          out[n][0][1] += product;
        }
        for (int c=i+j+1; c<=k; c++)
        {
          R product=factor;
          for (int beta=0; beta<j; beta++)
            product *= (x[1]-pos[beta])/(pos[j]-pos[beta]);
          for (int gamma=i+j+1; gamma<=k; gamma++)
            if (gamma==c)
              product *= -1.0/(pos[gamma]-pos[i]-pos[j]);
            else
              product *= (pos[gamma]-x[0]-x[1])/(pos[gamma]-pos[i]-pos[j]);
          out[n][0][1] += product;
        }

        n++;
      }

    //  for (int i=0; i<N; i++)
    //    std::cout << i << " " << out[i][0][0] << " " << out[i][0][1] << std::endl;
  }

  //! \brief Local interpolation of a function
  template<typename E, typename F, typename C>
  void interpolate (const E& e, const F& f, std::vector<C>& out) const
  {
    typename Traits::DomainType x;
    typename Traits::RangeType y;
    out.resize(N);
    int n=0;
    for (int j=0; j<=k; j++)
      for (int i=0; i<=k-j; i++)
      {
        x[0] = pos[i]; x[1] = pos[j];
        f.eval_local(e,x,y);
        out[n] = y;
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


//Specialization for k=0
template<class D, class R>
class Pk2DLocalBasis<D,R,0> :
  public C1LocalBasisInterface<
      C1LocalBasisTraits<D,2,Dune::FieldVector<D,2>,R,1,Dune::FieldVector<R,1>,
          Dune::FieldVector<Dune::FieldVector<R,2>,1> >,
      Pk2DLocalBasis<D,R,0> >
{
public:
  typedef C1LocalBasisTraits<D,2,Dune::FieldVector<D,2>,R,1,Dune::FieldVector<R,1>,
      Dune::FieldVector<Dune::FieldVector<R,2>,1> > Traits;

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
                    std::vector<typename Traits::JacobianType>& out) const          // return value
  {
    out.resize(1);
    out[0][0][0] = 0; out[0][0][1] = 0;
  }

  // local interpolation of a function
  template<typename E, typename F, typename C>
  void interpolate (const E& e, const F& f, std::vector<C>& out) const
  {
    typename Traits::DomainType x;
    typename Traits::RangeType y;
    x[0] = 1.0/3.0; x[1] = 1.0/3.0; f.eval_local(e,x,y); out[0] = y;
  }

  unsigned int order () const
  {
    return 0;
  }
};

#endif
