// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PK3DLOCALBASIS_HH
#define DUNE_PK3DLOCALBASIS_HH

#include "../common/localbasis.hh"

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
  class Pk3DLocalBasis :
    public C1LocalBasisInterface<
        C1LocalBasisTraits<D,3,Dune::FieldVector<D,3>,R,1,Dune::FieldVector<R,1>,
            Dune::FieldVector<Dune::FieldVector<R,3>,1> >
#ifndef DUNE_VIRTUAL_SHAPEFUNCTIONS
        ,Pk3DLocalBasis<D,R,k>
#endif
        >
  {
  public:
    enum {N = (k+1)*(k+2)*(k+3)/6};
    enum {O = k};

    typedef C1LocalBasisTraits<D,3,Dune::FieldVector<D,3>,R,1,Dune::FieldVector<R,1>,
        Dune::FieldVector<Dune::FieldVector<R,3>,1> > Traits;

    //! \brief Standard constructor
    Pk3DLocalBasis ()
    {
      for (unsigned int i=0; i<=k; i++)
        pos[i] = D(i)/k;
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
      unsigned int n = 0;
      for (unsigned int i3 = 0; i3 <= k; ++i3)
        for (unsigned int i2 = 0; i2 <= k - i3; ++i2)
          for (unsigned int i1 = 0; i1 <= k - i2 - i3; ++i1)
          {
            out[n] = 1.0;
            for (unsigned int alpha=0; alpha<i1; alpha++)
              out[n] *= (x[0]-pos[alpha])/(pos[i1]-pos[alpha]);
            for (unsigned int beta=0; beta<i2; beta++)
              out[n] *= (x[1]-pos[beta])/(pos[i2]-pos[beta]);
            for (unsigned int gamma=0; gamma<i3; gamma++)
              out[n] *= (x[2]-pos[gamma])/(pos[i3]-pos[gamma]);
            for (unsigned int delta=i1+i2+i3+1; delta<=k; delta++)
              out[n] *= (pos[delta]-x[0]-x[1]-x[2])/(pos[delta]-pos[i1]-pos[i2]-pos[i3]);
            n++;
          }
    }

    //! \brief Evaluate Jacobian of all shape functions
    inline void
    evaluateJacobian (const typename Traits::DomainType& x,         // position
                      std::vector<typename Traits::JacobianType>& out) const      // return value
    {
      DUNE_THROW(Dune::NotImplemented,
                 "Jacobian for P_k 3D shape functions not implemented yet");        out.resize(N);
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
  class Pk3DLocalBasis<D,R,0> :
    public C1LocalBasisInterface<
        C1LocalBasisTraits<D,3,Dune::FieldVector<D,3>,R,1,Dune::FieldVector<R,1>,
            Dune::FieldVector<Dune::FieldVector<R,3>,1> >
#ifndef DUNE_VIRTUAL_SHAPEFUNCTIONS
        ,Pk3DLocalBasis<D,R,0>
#endif
        >
  {
  public:
    typedef C1LocalBasisTraits<D,3,Dune::FieldVector<D,3>,R,1,Dune::FieldVector<R,1>,
        Dune::FieldVector<Dune::FieldVector<R,3>,1> > Traits;

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

    // local interpolation of a function
    template<typename E, typename F, typename C>
    void interpolate (const E& e, const F& f, std::vector<C>& out) const
    {
      typename Traits::DomainType x;
      typename Traits::RangeType y;
      x[0] = 1.0/3.0;
      x[1] = 1.0/3.0;
      x[2] = 1.0/3.0;
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
