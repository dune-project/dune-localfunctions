// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_RT0Q3DALL_HH
#define DUNE_RT0Q3DALL_HH

#include <cstddef>
#include <vector>

#include "../common/localbasis.hh"
#include "../common/localinterpolation.hh"
#include "../common/localcoefficients.hh"

namespace Dune
{
  /**@ingroup LocalBasisImplementation
         \brief Lowest order Raviart-Thomas shape functions on the reference hexahedron.

         - <tt>D</tt>: Type to represent the field in the domain.
         - <tt>R</tt>: Type to represent the field in the range.

         \nosubgrouping
   */
  template<class D, class R>
  class RT0Q3DLocalBasis
  {
  public:
    typedef C1LocalBasisTraits<D,3,Dune::FieldVector<D,3>,R,3,Dune::FieldVector<R,3>,
        Dune::FieldVector<Dune::FieldVector<R,3>,3> > Traits;

    //! \brief Standard constructor
    RT0Q3DLocalBasis ()
    {
      sign0 = sign1 = sign2 = sign3 = sign4 = sign5 = 1.0;
    }

    //! \brief Make set numer s, where 0<=s<64
    RT0Q3DLocalBasis (unsigned int s)
    {
      sign0 = sign1 = sign2 = sign3 = sign4 = sign5 = 1.0;
      if (s&1) sign0 = -1.0;
      if (s&2) sign1 = -1.0;
      if (s&4) sign2 = -1.0;
      if (s&8) sign3 = -1.0;
      if (s&16) sign4 = -1.0;
      if (s&32) sign5 = -1.0;
    }

    //! \brief number of shape functions
    unsigned int size () const
    {
      return 6;
    }

    //! \brief Evaluate all shape functions
    inline void evaluateFunction (const typename Traits::DomainType& in,
                                  std::vector<typename Traits::RangeType>& out) const
    {
      out.resize(6);
      out[0][0] = sign0*(in[0]-1.0); out[0][1]=0.0;               out[0][2]=0.0;
      out[1][0] = sign1*(in[0]);     out[1][1]=0.0;               out[1][2]=0.0;
      out[2][0] = 0.0;               out[2][1]=sign2*(in[1]-1.0); out[2][2]=0.0;
      out[3][0] = 0.0;               out[3][1]=sign3*(in[1]);     out[3][2]=0.0;
      out[4][0] = 0.0;               out[4][1]=0.0;               out[4][2]=sign4*(in[2]-1.0);
      out[5][0] = 0.0;               out[5][1]=0.0;               out[5][2]=sign5*(in[2]);
    }

    //! \brief Evaluate Jacobian of all shape functions
    inline void
    evaluateJacobian (const typename Traits::DomainType& in,             // position
                      std::vector<typename Traits::JacobianType>& out) const                          // return value
    {
      out.resize(6);
      out[0][0][0] = sign0;       out[0][0][1] = 0;      out[0][0][2] = 0;
      out[0][1][0] = 0;           out[0][1][1] = 0;      out[0][1][2] = 0;
      out[0][2][0] = 0;           out[0][2][1] = 0;      out[0][2][2] = 0;

      out[1][0][0] = sign1;       out[1][0][1] = 0;      out[1][0][2] = 0;
      out[1][1][0] = 0;           out[1][1][1] = 0;      out[1][1][2] = 0;
      out[1][2][0] = 0;           out[1][2][1] = 0;      out[1][2][2] = 0;

      out[2][0][0] = 0;           out[2][0][1] = 0;      out[2][0][2] = 0;
      out[2][1][0] = 0;           out[2][1][1] = sign2;  out[2][1][2] = 0;
      out[2][2][0] = 0;           out[2][2][1] = 0;      out[2][2][2] = 0;

      out[3][0][0] = 0;           out[3][0][1] = 0;      out[3][0][2] = 0;
      out[3][1][0] = 0;           out[3][1][1] = sign3;  out[3][1][2] = 0;
      out[3][2][0] = 0;           out[3][2][1] = 0;      out[3][2][2] = 0;

      out[4][0][0] = 0;           out[4][0][1] = 0;      out[4][0][2] = 0;
      out[4][1][0] = 0;           out[4][1][1] = 0;      out[4][1][2] = 0;
      out[4][2][0] = 0;           out[4][2][1] = 0;      out[4][2][2] = sign4;

      out[5][0][0] = 0;           out[5][0][1] = 0;      out[5][0][2] = 0;
      out[5][1][0] = 0;           out[5][1][1] = 0;      out[5][1][2] = 0;
      out[5][2][0] = 0;           out[5][2][1] = 0;      out[5][2][2] = sign5;
    }

    //! \brief Polynomial order of the shape functions
    unsigned int order () const
    {
      return 1;
    }

  private:
    R sign0, sign1, sign2, sign3, sign4, sign5;
  };


  /**@ingroup LocalInterpolationImplementation
         \brief Lowest order Raviart-Thomas shape functions on the reference hexahedron.

         - <tt>LB</tt>: corresponding LocalBasis giving traits

         \nosubgrouping
   */
  template<class LB>
  class RT0Q3DLocalInterpolation
  {
  public:

    //! \brief Standard constructor
    RT0Q3DLocalInterpolation ()
    {}

    //! \brief Make set numer s, where 0<=s<64
    RT0Q3DLocalInterpolation (unsigned int s)
    {
      sign0 = sign1 = sign2 = sign3 = sign4 = sign5 = 1.0;
      if (s&1) sign0 *= -1.0;
      if (s&2) sign1 *= -1.0;
      if (s&4) sign2 *= -1.0;
      if (s&8) sign3 *= -1.0;
      if (s&16) sign4 *= -1.0;
      if (s&32) sign5 *= -1.0;

      m0[0] = 0.0; m0[1] = 0.5; m0[2] = 0.5;
      m1[0] = 1.0; m1[1] = 0.5; m1[2] = 0.5;
      m2[0] = 0.5; m2[1] = 0.0; m2[2] = 0.5;
      m3[0] = 0.5; m3[1] = 1.0; m3[2] = 0.5;
      m4[0] = 0.5; m4[1] = 0.5; m4[2] = 0.0;
      m5[0] = 0.5; m5[1] = 0.5; m5[2] = 1.0;

      n0[0] = -1.0; n0[1] =  0.0; n0[2] = 0.0;
      n1[0] =  1.0; n1[1] =  0.0; n1[2] = 0.0;
      n2[0] =  0.0; n2[1] = -1.0; n2[2] = 0.0;
      n3[0] =  0.0; n3[1] =  1.0; n3[2] = 0.0;
      n4[0] =  0.0; n4[1] =  0.0; n4[2] =-1.0;
      n5[0] =  0.0; n5[1] =  0.0; n5[2] = 1.0;
    }

    template<typename F, typename C>
    void interpolate (const F& f, std::vector<C>& out) const
    {
      // f gives v*outer normal at a point on the edge!
      typename F::Traits::RangeType y;

      out.resize(6);

      f.evaluate(m0,y); out[0] = (y[0]*n0[0]+y[1]*n0[1]+y[2]*n0[2])*sign0;
      f.evaluate(m1,y); out[1] = (y[0]*n1[0]+y[1]*n1[1]+y[2]*n1[2])*sign1;
      f.evaluate(m2,y); out[2] = (y[0]*n2[0]+y[1]*n2[1]+y[2]*n2[2])*sign2;
      f.evaluate(m3,y); out[3] = (y[0]*n3[0]+y[1]*n3[1]+y[2]*n3[2])*sign3;
      f.evaluate(m4,y); out[4] = (y[0]*n4[0]+y[1]*n4[1]+y[2]*n4[2])*sign4;
      f.evaluate(m5,y); out[5] = (y[0]*n5[0]+y[1]*n5[1]+y[2]*n5[2])*sign5;
    }

  private:
    typename LB::Traits::RangeFieldType sign0,sign1,sign2,sign3,sign4,sign5;
    typename LB::Traits::DomainType m0,m1,m2,m3,m4,m5;
    typename LB::Traits::DomainType n0,n1,n2,n3,n4,n5;
  };

  /**@ingroup LocalLayoutImplementation
         \brief Layout map for RT0 elements on quadrilaterals

         \nosubgrouping
     \implements Dune::LocalCoefficientsVirtualImp
   */
  class RT0Q3DLocalCoefficients
  {
  public:
    //! \brief Standard constructor
    RT0Q3DLocalCoefficients () : li(6)
    {
      for (std::size_t i=0; i<6; i++)
        li[i] = LocalKey(i,1,0);
    }

    //! number of coefficients
    std::size_t size () const
    {
      return 6;
    }

    //! get i'th index
    const LocalKey& localKey (std::size_t i) const
    {
      return li[i];
    }

  private:
    std::vector<LocalKey> li;
  };

}
#endif
