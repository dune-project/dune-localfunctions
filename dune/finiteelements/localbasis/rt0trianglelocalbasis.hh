// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_RT0TRIANGLELOCALBASIS_HH
#define DUNE_RT0TRIANGLELOCALBASIS_HH

#include "localbasis.hh"

/**@ingroup LocalBasisImplementation
   \brief Lowest order Raviart-Thomas shape functions on the reference triangle.

   - <tt>D</tt>: Type to represent the field in the domain.
   - <tt>R</tt>: Type to represent the field in the range.

   \nosubgrouping
 */
template<class D, class R>
class RT0TriangleLocalBasis :
  public C1LocalBasisInterface<
      C1LocalBasisTraits<D,2,Dune::FieldVector<D,2>,R,2,Dune::FieldVector<R,2>,
          Dune::FieldVector<Dune::FieldVector<R,2>,2> >,
      RT0TriangleLocalBasis<D,R> >
{
public:
  typedef C1LocalBasisTraits<D,2,Dune::FieldVector<D,2>,R,2,Dune::FieldVector<R,2> > Traits;

  //! \brief Standard constructor
  RT0TriangleLocalBasis ()
  {
    sqrt2 = sqrt(2.0);
    sign0 = sign1 = sign2 = 1;
  }

  //! \brief Make set numer s, where 0<=s<8
  RT0TriangleLocalBasis (unsigned int s)
  {
    sqrt2 = sqrt(2.0);
    sign0 = sign1 = sign2 = 1;
    if (s&1==0) sign0 *= -1;
    if (s&2==0) sign1 *= -1;
    if (s&4==0) sign2 *= -1;
  }

  //! \brief number of shape functions
  unsigned int size () const
  {
    return 3;
  }

  //! \brief Evaluate all shape functions
  inline void evaluateFunction (const typename Traits::DomainType& in,
                                std::vector<typename Traits::RangeType>& out) const
  {
    out.resize(3);
    out[0][0] = sign0*in[0]*sqrt2; out[0][1]=sign0*in[1]*sqrt2;
    out[1][0] = sign1*(in[0]-1);   out[1][1]=sign1*(in[1]);
    out[2][0] = sign2*in[0];       out[2][1]=sign2*in[1]-1;
  }

  //! \brief Evaluate Jacobian of all shape functions
  inline void
  evaluateJacobian (const typename Traits::DomainType& in,         // position
                    std::vector<typename Traits::JacobianType>& out) const          // return value
  {
    out.resize(3);
    out[0][0][0] = sign0*sqrt2; out[0][1][0] = 0; out[0][0][1] = 0; out[0][1][1] = sign0*sqrt2;
    out[1][0][0] = 1;           out[1][1][0] = 0; out[1][0][1] = 0; out[1][1][1] = 1;
    out[2][0][0] = 1;           out[2][1][0] = 0; out[2][0][1] = 0; out[2][1][1] = 1;
  }

  //! \brief Local interpolation of a function
  template<typename E, typename F, typename C>
  void interpolate (const E& e, const F& f, std::vector<C>& out) const
  {
    typename Traits::DomainType x;
    typename Traits::RangeType y;

    out.resize(3);
    DUNE_THROW(Dune::Exception,"interpolate for RT0 not implemented");
  }

  //! \brief Polynomial order of the shape functions
  unsigned int order () const
  {
    return 1;
  }

private:
  R sqrt2;
  R sign0, sign1, sign2;
};

#endif
