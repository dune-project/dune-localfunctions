// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS3_CUBE2D_LOCALBASIS_HH
#define DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS3_CUBE2D_LOCALBASIS_HH

#include <vector>

#include <dune/common/fmatrix.hh>

#include "../../common/localbasis.hh"

namespace Dune
{
  /**
   * \ingroup LocalBasisImplementation
   * \brief Second order Raviart-Thomas shape functions on the reference quadrilateral.
   *
   * \tparam D Type to represent the field in the domain.
   * \tparam R Type to represent the field in the range.
   *
   * \nosubgrouping
   */
  template<class D, class R>
  class RT3Cube2DLocalBasis
  {

  public:
    typedef LocalBasisTraits<D,2,Dune::FieldVector<D,2>,R,2,Dune::FieldVector<R,2>,
        Dune::FieldMatrix<R,2,2>, 0 > Traits;

    //! \brief Standard constructor
    RT3Cube2DLocalBasis ()
    {
      sign0 = sign1 = sign2 = sign3 = 1.0;
    }

    /**
     * \brief Make set number s, where 0 <= s < 16
     *
     * \param s Edge orientation indicator
     */
    RT3Cube2DLocalBasis (unsigned int s)
    {
      sign0 = sign1 = sign2 = sign3 = 1.0;
      if (s & 1)
      {
        sign0 = -1.0;
      }
      if (s & 2)
      {
        sign1 = -1.0;
      }
      if (s & 4)
      {
        sign2 = -1.0;
      }
      if (s & 8)
      {
        sign3 = -1.0;
      }
    }

    //! \brief number of shape functions
    unsigned int size () const
    {
      return 40;
    }

    /**
     * \brief Evaluate all shape functions
     *
     * \param in Position
     * \param out return value
     */
    inline void evaluateFunction (const typename Traits::DomainType& in,
                                  std::vector<typename Traits::RangeType>& out) const
    {
      out.resize(40);
      double tmp1=-1.0+16.0*in[0]-60.0*pow(in[0],2)+80.0*pow(in[0],3)-35.0*pow(in[0],4);
      double tmp2=(1.0-16.0*in[0]+60.0*pow(in[0],2)-80.0*pow(in[0],3)+35.0*pow(in[0],4));
      double tmp3=(-1.0+2.0*in[1]);
      double tmp4=(1.0-6.0*in[1]+6.0*pow(in[1],2));
      double tmp5=(-1.0+12.0*in[1]-30.0*pow(in[1],2)+20.0*pow(in[1],3));
      double tmp6=in[0]*(-4.0+30.0*in[0]-60.0*pow(in[0],2)+35.0*pow(in[0],3));
      double tmp7=-1.0+16.0*in[1]-60.0*pow(in[1],2)+80.0*pow(in[1],3)-35.0*pow(in[1],4);
      double tmp8=(1.0-16.0*in[1]+60.0*pow(in[1],2)-80.0*pow(in[1],3)+35.0*pow(in[1],4));
      double tmp9=(-1.0+2.0*in[0]);
      double tmp10=(1.0-6.0*in[0]+6.0*pow(in[0],2));
      double tmp11=(-1.0+12.0*in[0]-30.0*pow(in[0],2)+20.0*pow(in[0],3));
      double tmp12=in[1]*(-4.0+30.0*in[1]-60.0*pow(in[1],2)+35.0*pow(in[1],3));
      double tmp13=in[0]*(2.0-9.0*in[0]+14.0*pow(in[0],2)-7.0*pow(in[0],3));
      double tmp14=in[0]*(-2.0+9.0*in[0]-14.0*pow(in[0],2)+7.0*pow(in[0],3));
      double tmp15=in[0]*(1.0-3.0*in[0]+2.0*pow(in[0],2));
      double tmp16=in[0]*(-1.0+6.0*in[0]-10.0*pow(in[0],2)+5.0*pow(in[0],3));
      double tmp17=in[1]*(2.0-9.0*in[1]+14.0*pow(in[1],2)-7.0*pow(in[1],3));
      double tmp18=in[1]*(1.0-3.0*in[1]+2.0*pow(in[1],2));
      double tmp19=in[1]*(-1.0+6.0*in[1]-10.0*pow(in[1],2)+5.0*pow(in[1],3));
      double tmp20=in[1]*(-2.0+9.0*in[1]-14.0*pow(in[1],2)+7.0*pow(in[1],3));

      out[0][0]=sign0*tmp1;
      out[0][1]=0;
      out[1][0]=(-3.0*tmp2*tmp3);
      out[1][1]=0;
      out[2][0]=sign0*(-5.0*tmp2*tmp4);
      out[2][1]=0;
      out[3][0]=(-7.0*tmp2*tmp5);
      out[3][1]=0;

      out[4][0]=sign1*tmp6;
      out[4][1]=0;
      out[5][0]=(-3.0*tmp6*tmp3);
      out[5][1]=0;
      out[6][0]=sign1*(5.0*tmp6*tmp4);
      out[6][1]=0;
      out[7][0]=(-7.0*tmp6*tmp5);
      out[7][1]=0;

      out[8][0]=0;
      out[8][1]=sign2*tmp7;
      out[9][0]=0;
      out[9][1]=3.0*tmp9*tmp8;
      out[10][0]=0;
      out[10][1]=sign2*(-5.0*tmp10*tmp8);
      out[11][0]=0;
      out[11][1]=7.0*tmp11*tmp8;

      out[12][0]=0;
      out[12][1]=sign3*tmp12;
      out[13][0]=0;
      out[13][1]=3.0*tmp9*tmp12;
      out[14][0]=0;
      out[14][1]=sign3*5.0*tmp10*tmp12;
      out[15][0]=0;
      out[15][1]=7.0*tmp11*tmp12;

      out[16][0]=10.0*tmp13;
      out[16][1]=0;
      out[17][0]=-30.0*tmp14*tmp3;
      out[17][1]=0;
      out[18][0]=-50.0*tmp14*tmp4;
      out[18][1]=0;
      out[19][0]=-70.0*tmp14*tmp5;
      out[19][1]=0;
      out[20][0]=-30.0*tmp15;
      out[20][1]=0;
      out[21][0]=-90.0*tmp15*tmp3;
      out[21][1]=0;
      out[22][0]=-150.0*tmp15*tmp4;
      out[22][1]=0;
      out[23][0]=-210.0*tmp15*tmp5;
      out[23][1]=0;
      out[24][0]=-70.0*tmp16;
      out[24][1]=0;
      out[25][0]=-210.0*tmp16*tmp3;
      out[25][1]=0;
      out[26][0]=-350.0*tmp16*tmp4;
      out[26][1]=0;
      out[27][0]=-490.0*tmp16*tmp5;
      out[27][1]=0;
      out[28][0]=0;
      out[28][1]=10.0*tmp17;
      out[29][0]=0;
      out[29][1]=-30.0*tmp18;
      out[30][0]=0;
      out[30][1]=-70.0*tmp19;
      out[31][0]=0;
      out[31][1]=-30.0*tmp9*tmp20;
      out[32][0]=0;
      out[32][1]=-90.0*tmp9*tmp18;
      out[33][0]=0;
      out[33][1]=-210.0*tmp9*tmp19;
      out[34][0]=0;
      out[34][1]=-50.0*tmp10*tmp20;
      out[35][0]=0;
      out[35][1]=-150.0*tmp10*tmp18;
      out[36][0]=0;
      out[36][1]=-350.0*tmp10*tmp19;
      out[37][0]=0;
      out[37][1]=-70.0*tmp11*tmp20;
      out[38][0]=0;
      out[38][1]=-210.0*tmp11*tmp18;
      out[39][0]=0;
      out[39][1]=-490.0*tmp11*tmp19;
    }

    /**
     * \brief Evaluate Jacobian of all shape functions
     *
     * \param in Position
     * \param out return value
     */
    inline void evaluateJacobian (const typename Traits::DomainType& /*in*/,
                                  std::vector<typename Traits::JacobianType>& out) const
    {
      out.resize(40);
      DUNE_THROW(NotImplemented, "Jacobian is not implemented");
    }

    //! \brief Evaluate partial derivatives of all shape functions
    inline void partial (const std::array<unsigned int, 2>& order,
                         const typename Traits::DomainType& in,         // position
                         std::vector<typename Traits::RangeType>& out) const      // return value
    {
      auto totalOrder = std::accumulate(order.begin(), order.end(), 0);
      if (totalOrder == 0) {
        evaluateFunction(in, out);
      } else {
        DUNE_THROW(NotImplemented, "Desired derivative order is not implemented");
      }
    }

    //! \brief Evaluate partial derivatives of all shape functions, \deprecated
    template <std::size_t dOrder>
    inline void evaluate (const std::array<int, dOrder>& directions,
                          const typename Traits::DomainType& in,         // position
                          std::vector<typename Traits::RangeType>& out) const      // return value
    {
      std::array<unsigned int, 2> order;
      Impl::directions2order(directions, order);
      partial(order, in, out);
    }

    //! \brief Polynomial order of the shape functions
    unsigned int order () const
    {
      return 7;
    }

  private:
    R sign0, sign1, sign2, sign3;
  };
}

#endif // DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS3_CUBE2D_LOCALBASIS_HH
