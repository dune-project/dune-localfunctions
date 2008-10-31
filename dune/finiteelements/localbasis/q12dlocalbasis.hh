// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_Q12DLOCALBASIS_HH
#define DUNE_Q12DLOCALBASIS_HH

#include "localbasis.hh"

/**@ingroup LocalBasisImplementation
   \brief Lagrange shape functions of order 1 on the reference quadrilateral.

   Also known as \f$Q^1\f$.

   - <tt>D</tt>: Type to represent the field in the domain.
   - <tt>R</tt>: Type to represent the field in the range.

   \nosubgrouping
 */
template<class D, class R>
class Q12DLocalBasis :
  public C1LocalBasisInterface<
      C1LocalBasisTraits<D,2,Dune::FieldVector<D,2>,R,1,Dune::FieldVector<R,1>,
          Dune::FieldVector<Dune::FieldVector<R,2>,1> >,
      Q12DLocalBasis<D,R> >
{
public:
  typedef C1LocalBasisTraits<D,2,Dune::FieldVector<D,2>,R,1,Dune::FieldVector<R,1>,
      Dune::FieldVector<Dune::FieldVector<R,2>,1> > Traits;

  //! \brief number of shape functions
  unsigned int size () const
  {
    return 4;
  }

  //! \brief Evaluate all shape functions
  inline void evaluateFunction (const typename Traits::DomainType& in,
                                std::vector<typename Traits::RangeType>& out) const
  {
    out.resize(4);
    out[0] = (1-in[0])*(1-in[1]);
    out[1] = (  in[0])*(1-in[1]);
    out[2] = (1-in[0])*(  in[1]);
    out[3] = (  in[0])*(  in[1]);
  }

  //! \brief Evaluate Jacobian of all shape functions
  inline void
  evaluateJacobian (const typename Traits::DomainType& in,         // position
                    std::vector<typename Traits::JacobianType>& out) const          // return value
  {
    out.resize(4);
    out[0][0][0] = in[1]-1; out[0][0][1] = in[0]-1;
    out[1][0][0] = 1-in[1]; out[1][0][1] = -in[0];
    out[2][0][0] =  -in[1]; out[2][0][1] = 1-in[0];
    out[3][0][0] =   in[1]; out[3][0][1] = in[0];
  }

  //! \brief Local interpolation of a function
  template<typename E, typename F, typename C>
  void interpolate (const E& e, const F& f, std::vector<C>& out) const
  {
    typename Traits::DomainType x;
    typename Traits::RangeType y;

    out.resize(4);
    x[0] = 0.0; x[1] = 0.0; f.eval_local(e,x,y); out[0] = y;
    x[0] = 1.0; x[1] = 0.0; f.eval_local(e,x,y); out[1] = y;
    x[0] = 0.0; x[1] = 1.0; f.eval_local(e,x,y); out[2] = y;
    x[0] = 1.0; x[1] = 1.0; f.eval_local(e,x,y); out[3] = y;
  }

  //! \brief Polynomial order of the shape functions
  unsigned int order () const
  {
    return 1;
  }
};

#endif
