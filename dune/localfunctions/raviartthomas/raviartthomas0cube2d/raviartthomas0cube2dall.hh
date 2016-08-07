// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS0_CUBE2D_ALL_HH
#define DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS0_CUBE2D_ALL_HH

#include <cstddef>
#include <numeric>
#include <vector>

#include <dune/common/fmatrix.hh>

#include <dune/localfunctions/common/localbasis.hh>
#include <dune/localfunctions/common/localkey.hh>

namespace Dune
{
  /**@ingroup LocalBasisImplementation
         \brief Lowest order Raviart-Thomas shape functions on the reference quadrilateral.

         \tparam D Type to represent the field in the domain.
         \tparam R Type to represent the field in the range.

         \nosubgrouping
   */
  template<class D, class R>
  class RT0Cube2DLocalBasis
  {
  public:
    typedef LocalBasisTraits<D,2,Dune::FieldVector<D,2>,R,2,Dune::FieldVector<R,2>,
        Dune::FieldMatrix<R,2,2> > Traits;

    //! \brief Standard constructor
    RT0Cube2DLocalBasis ()
    {
      sign0 = sign1 = sign2 = sign3 = 1.0;
    }

    //! \brief Make set numer s, where 0<=s<16
    RT0Cube2DLocalBasis (unsigned int s)
    {
      sign0 = sign1 = sign2 = sign3 = 1.0;
      if (s&1) sign0 = -1.0;
      if (s&2) sign1 = -1.0;
      if (s&4) sign2 = -1.0;
      if (s&8) sign3 = -1.0;
    }

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
      out[0][0] = sign0*(in[0]-1.0); out[0][1]=0.0;
      out[1][0] = sign1*(in[0]);     out[1][1]=0.0;
      out[2][0] = 0.0;               out[2][1]=sign2*(in[1]-1.0);
      out[3][0] = 0.0;               out[3][1]=sign3*(in[1]);
    }

    //! \brief Evaluate Jacobian of all shape functions
    inline void
    evaluateJacobian (const typename Traits::DomainType& in,             // position
                      std::vector<typename Traits::JacobianType>& out) const                          // return value
    {
      out.resize(4);
      out[0][0][0] = sign0;       out[0][0][1] = 0;
      out[0][1][0] = 0;           out[0][1][1] = 0;

      out[1][0][0] = sign1;       out[1][0][1] = 0;
      out[1][1][0] = 0;           out[1][1][1] = 0;

      out[2][0][0] = 0;           out[2][0][1] = 0;
      out[2][1][0] = 0;           out[2][1][1] = sign2;

      out[3][0][0] = 0;           out[3][0][1] = 0;
      out[3][1][0] = 0;           out[3][1][1] = sign3;
    }

    //! \brief Evaluate partial derivatives of all shape functions
    void partial (const std::array<unsigned int, 2>& order,
                  const typename Traits::DomainType& in,         // position
                  std::vector<typename Traits::RangeType>& out) const      // return value
    {
      auto totalOrder = std::accumulate(order.begin(), order.end(), 0);
      if (totalOrder == 0) {
        evaluateFunction(in, out);
      } else if (totalOrder == 1) {
        auto const direction = std::distance(order.begin(), std::find(order.begin(), order.end(), 1));
        out.resize(size());

        for (std::size_t i = 0; i < size(); ++i)
          out[i][0] = out[i][1] = 0;

        switch (direction) {
        case 0:
          out[0][0] = sign0;
          out[1][0] = sign1;
          break;
        case 1:
          out[2][1] = sign2;
          out[3][1] = sign3;
          break;
        default:
          DUNE_THROW(RangeError, "Component out of range.");
        }
      } else {
        out.resize(size());
        for (std::size_t i = 0; i < size(); ++i)
          for (std::size_t j = 0; j < 2; ++j)
            out[i][j] = 0;
      }

    }

    //! \brief Polynomial order of the shape functions
    unsigned int order () const
    {
      return 1;
    }

  private:
    R sign0, sign1, sign2, sign3;
  };


  /**@ingroup LocalInterpolationImplementation
         \brief Lowest order Raviart-Thomas shape functions on the reference quadrilateral.

         \tparam LB corresponding LocalBasis giving traits

         \nosubgrouping
   */
  template<class LB>
  class RT0Cube2DLocalInterpolation
  {
  public:

    //! \brief Standard constructor
    RT0Cube2DLocalInterpolation ()
    {
      sign0 = sign1 = sign2 = sign3 = 1.0;
    }

    //! \brief Make set numer s, where 0<=s<8
    RT0Cube2DLocalInterpolation (unsigned int s)
    {
      sign0 = sign1 = sign2 = sign3 = 1.0;
      if (s&1) sign0 *= -1.0;
      if (s&2) sign1 *= -1.0;
      if (s&4) sign2 *= -1.0;
      if (s&8) sign3 *= -1.0;

      m0[0] = 0.0; m0[1] = 0.5;
      m1[0] = 1.0; m1[1] = 0.5;
      m2[0] = 0.5; m2[1] = 0.0;
      m3[0] = 0.5; m3[1] = 1.0;

      n0[0] = -1.0; n0[1] =  0.0;
      n1[0] =  1.0; n1[1] =  0.0;
      n2[0] =  0.0; n2[1] = -1.0;
      n3[0] =  0.0; n3[1] =  1.0;
    }

    template<typename F, typename C>
    void interpolate (const F& f, std::vector<C>& out) const
    {
      // f gives v*outer normal at a point on the edge!
      typename F::Traits::RangeType y;

      out.resize(4);

      f.evaluate(m0,y); out[0] = (y[0]*n0[0]+y[1]*n0[1])*sign0;
      f.evaluate(m1,y); out[1] = (y[0]*n1[0]+y[1]*n1[1])*sign1;
      f.evaluate(m2,y); out[2] = (y[0]*n2[0]+y[1]*n2[1])*sign2;
      f.evaluate(m3,y); out[3] = (y[0]*n3[0]+y[1]*n3[1])*sign3;
    }

  private:
    typename LB::Traits::RangeFieldType sign0,sign1,sign2,sign3;
    typename LB::Traits::DomainType m0,m1,m2,m3;
    typename LB::Traits::DomainType n0,n1,n2,n3;
  };

  /**@ingroup LocalLayoutImplementation
         \brief Layout map for RT0 elements on quadrilaterals

         \nosubgrouping
     \implements Dune::LocalCoefficientsVirtualImp
   */
  class RT0Cube2DLocalCoefficients
  {
  public:
    //! \brief Standard constructor
    RT0Cube2DLocalCoefficients () : li(4)
    {
      for (std::size_t i=0; i<4; i++)
        li[i] = LocalKey(i,1,0);
    }

    //! number of coefficients
    std::size_t size () const
    {
      return 4;
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
#endif // DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS0_CUBE2D_ALL_HH
