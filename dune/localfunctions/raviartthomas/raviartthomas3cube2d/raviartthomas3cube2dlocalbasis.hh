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
        Dune::FieldMatrix<R,2,2> > Traits;

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
      auto const& x = in[0], y = in[1];

      auto tmp1 = - x*(x*(x*(35*x - 80) + 60) - 16) - 1;
      auto tmp2 = x*(x*(x*(35*x - 80) + 60) - 16) + 1;
      auto tmp3 = 2*y - 1;
      auto tmp4 = y*(6*y - 6) + 1;
      auto tmp5 = y*(y*(20*y - 30) + 12) - 1;
      auto tmp6 = x*(x*(x*(35*x - 60) + 30) - 4);
      auto tmp7 = - y*(y*(y*(35*y - 80) + 60) - 16) - 1;
      auto tmp8 = y*(y*(y*(35*y - 80) + 60) - 16) + 1;
      auto tmp9 = 2*x - 1;
      auto tmp10 = x*(6*x - 6) + 1;
      auto tmp11 = x*(x*(20*x - 30) + 12) - 1;
      auto tmp12 = y*(y*(y*(35*y - 60) + 30) - 4);
      auto tmp13 = -x*(x*(x*(7*x - 14) + 9) - 2);
      auto tmp14 = x*(x*(x*(7*x - 14) + 9) - 2);
      auto tmp15 = x*(x*(2*x - 3) + 1);
      auto tmp16 = x*(x*(x*(5*x - 10) + 6) - 1);
      auto tmp17 = -y*(y*(y*(7*y - 14) + 9) - 2);
      auto tmp18 = y*(y*(2*y - 3) + 1);
      auto tmp19 = y*(y*(y*(5*y - 10) + 6) - 1);
      auto tmp20 = y*(y*(y*(7*y - 14) + 9) - 2);

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
    inline void evaluateJacobian (const typename Traits::DomainType& in,
                                  std::vector<typename Traits::JacobianType>& out) const
    {
      out.resize(40);
      auto const& x = in[0], y = in[1];

      auto tmp1 = - x*(x*(x*(35*x - 80) + 60) - 16) - 1;
      auto tmp2 = x*(x*(x*(35*x - 80) + 60) - 16) + 1;
      auto tmp3 = 2*y - 1;
      auto tmp4 = y*(6*y - 6) + 1;
      auto tmp5 = y*(y*(20*y - 30) + 12) - 1;
      auto tmp6 = x*(x*(x*(35*x - 60) + 30) - 4);
      auto tmp7 = - y*(y*(y*(35*y - 80) + 60) - 16) - 1;
      auto tmp8 = y*(y*(y*(35*y - 80) + 60) - 16) + 1;
      auto tmp9 = 2*x - 1;
      auto tmp10 = x*(6*x - 6) + 1;
      auto tmp11 = x*(x*(20*x - 30) + 12) - 1;
      auto tmp12 = y*(y*(y*(35*y - 60) + 30) - 4);
      auto tmp13 = -x*(x*(x*(7*x - 14) + 9) - 2);
      auto tmp14 = x*(x*(x*(7*x - 14) + 9) - 2);
      auto tmp15 = x*(x*(2*x - 3) + 1);
      auto tmp16 = x*(x*(x*(5*x - 10) + 6) - 1);
      auto tmp17 = -y*(y*(y*(7*y - 14) + 9) - 2);
      auto tmp18 = y*(y*(2*y - 3) + 1);
      auto tmp19 = y*(y*(y*(5*y - 10) + 6) - 1);
      auto tmp20 = y*(y*(y*(7*y - 14) + 9) - 2);

      auto dxtmp1 = 16 - x*(x*(140*x - 240) + 120);
      auto dxtmp2 = x*(x*(140*x - 240) + 120) - 16;
      auto dytmp3 = 2;
      auto dytmp4 = 12*y - 6;
      auto dytmp5 = y*(60*y - 60) + 12;
      auto dxtmp6 = x*(x*(140*x - 180) + 60) - 4;
      auto dytmp7 = 16 - y*(y*(140*y - 240) + 120);
      auto dytmp8 = y*(y*(140*y - 240) + 120) - 16;
      auto dxtmp9 = 2;
      auto dxtmp10 = 12*x - 6;
      auto dxtmp11 = x*(60*x - 60) + 12;
      auto dytmp12 = y*(y*(140*y - 180) + 60) - 4;
      auto dxtmp13 = 2 - x*(x*(28*x - 42) + 18);
      auto dxtmp14 = x*(x*(28*x - 42) + 18) - 2;
      auto dxtmp15 = x*(6*x - 6) + 1;
      auto dxtmp16 = x*(x*(20*x - 30) + 12) - 1;
      auto dytmp17 = 2 - y*(y*(28*y - 42) + 18);
      auto dytmp18 = y*(6*y - 6) + 1;
      auto dytmp19 = y*(y*(20*y - 30) + 12) - 1;
      auto dytmp20 = y*(y*(28*y - 42) + 18) - 2;


      // x-component
      out[0][0][0]=sign0*dxtmp1;
      out[0][1][0]=0;
      out[1][0][0]=(-3.0*dxtmp2*tmp3);
      out[1][1][0]=0;
      out[2][0][0]=sign0*(-5.0*dxtmp2*tmp4);
      out[2][1][0]=0;
      out[3][0][0]=(-7.0*dxtmp2*tmp5);
      out[3][1][0]=0;

      out[4][0][0]=sign1*dxtmp6;
      out[4][1][0]=0;
      out[5][0][0]=(-3.0*dxtmp6*tmp3);
      out[5][1][0]=0;
      out[6][0][0]=sign1*(5.0*dxtmp6*tmp4);
      out[6][1][0]=0;
      out[7][0][0]=(-7.0*dxtmp6*tmp5);
      out[7][1][0]=0;

      out[8][0][0]=0;
      out[8][1][0]=0;
      out[9][0][0]=0;
      out[9][1][0]=3.0*dxtmp9*tmp8;
      out[10][0][0]=0;
      out[10][1][0]=sign2*(-5.0*dxtmp10*tmp8);
      out[11][0][0]=0;
      out[11][1][0]=7.0*dxtmp11*tmp8;

      out[12][0][0]=0;
      out[12][1][0]=0;
      out[13][0][0]=0;
      out[13][1][0]=3.0*dxtmp9*tmp12;
      out[14][0][0]=0;
      out[14][1][0]=sign3*5.0*dxtmp10*tmp12;
      out[15][0][0]=0;
      out[15][1][0]=7.0*dxtmp11*tmp12;

      out[16][0][0]=10.0*dxtmp13;
      out[16][1][0]=0;
      out[17][0][0]=-30.0*dxtmp14*tmp3;
      out[17][1][0]=0;
      out[18][0][0]=-50.0*dxtmp14*tmp4;
      out[18][1][0]=0;
      out[19][0][0]=-70.0*dxtmp14*tmp5;
      out[19][1][0]=0;
      out[20][0][0]=-30.0*dxtmp15;
      out[20][1][0]=0;
      out[21][0][0]=-90.0*dxtmp15*tmp3;
      out[21][1][0]=0;
      out[22][0][0]=-150.0*dxtmp15*tmp4;
      out[22][1][0]=0;
      out[23][0][0]=-210.0*dxtmp15*tmp5;
      out[23][1][0]=0;
      out[24][0][0]=-70.0*dxtmp16;
      out[24][1][0]=0;
      out[25][0][0]=-210.0*dxtmp16*tmp3;
      out[25][1][0]=0;
      out[26][0][0]=-350.0*dxtmp16*tmp4;
      out[26][1][0]=0;
      out[27][0][0]=-490.0*dxtmp16*tmp5;
      out[27][1][0]=0;
      out[28][0][0]=0;
      out[28][1][0]=0;
      out[29][0][0]=0;
      out[29][1][0]=0;
      out[30][0][0]=0;
      out[30][1][0]=0;
      out[31][0][0]=0;
      out[31][1][0]=-30.0*dxtmp9*tmp20;
      out[32][0][0]=0;
      out[32][1][0]=-90.0*dxtmp9*tmp18;
      out[33][0][0]=0;
      out[33][1][0]=-210.0*dxtmp9*tmp19;
      out[34][0][0]=0;
      out[34][1][0]=-50.0*dxtmp10*tmp20;
      out[35][0][0]=0;
      out[35][1][0]=-150.0*dxtmp10*tmp18;
      out[36][0][0]=0;
      out[36][1][0]=-350.0*dxtmp10*tmp19;
      out[37][0][0]=0;
      out[37][1][0]=-70.0*dxtmp11*tmp20;
      out[38][0][0]=0;
      out[38][1][0]=-210.0*dxtmp11*tmp18;
      out[39][0][0]=0;
      out[39][1][0]=-490.0*dxtmp11*tmp19;


      // y-component
      out[0][0][1]=0;
      out[0][1][1]=0;
      out[1][0][1]=(-3.0*tmp2*dytmp3);
      out[1][1][1]=0;
      out[2][0][1]=sign0*(-5.0*tmp2*dytmp4);
      out[2][1][1]=0;
      out[3][0][1]=(-7.0*tmp2*dytmp5);
      out[3][1][1]=0;

      out[4][0][1]=0;
      out[4][1][1]=0;
      out[5][0][1]=(-3.0*tmp6*dytmp3);
      out[5][1][1]=0;
      out[6][0][1]=sign1*(5.0*tmp6*dytmp4);
      out[6][1][1]=0;
      out[7][0][1]=(-7.0*tmp6*dytmp5);
      out[7][1][1]=0;

      out[8][0][1]=0;
      out[8][1][1]=sign2*dytmp7;
      out[9][0][1]=0;
      out[9][1][1]=3.0*tmp9*dytmp8;
      out[10][0][1]=0;
      out[10][1][1]=sign2*(-5.0*tmp10*dytmp8);
      out[11][0][1]=0;
      out[11][1][1]=7.0*tmp11*dytmp8;

      out[12][0][1]=0;
      out[12][1][1]=sign3*dytmp12;
      out[13][0][1]=0;
      out[13][1][1]=3.0*tmp9*dytmp12;
      out[14][0][1]=0;
      out[14][1][1]=sign3*5.0*tmp10*dytmp12;
      out[15][0][1]=0;
      out[15][1][1]=7.0*tmp11*dytmp12;

      out[16][0][1]=0;
      out[16][1][1]=0;
      out[17][0][1]=-30.0*tmp14*dytmp3;
      out[17][1][1]=0;
      out[18][0][1]=-50.0*tmp14*dytmp4;
      out[18][1][1]=0;
      out[19][0][1]=-70.0*tmp14*dytmp5;
      out[19][1][1]=0;
      out[20][0][1]=0;
      out[20][1][1]=0;
      out[21][0][1]=-90.0*tmp15*dytmp3;
      out[21][1][1]=0;
      out[22][0][1]=-150.0*tmp15*dytmp4;
      out[22][1][1]=0;
      out[23][0][1]=-210.0*tmp15*dytmp5;
      out[23][1][1]=0;
      out[24][0][1]=0;
      out[24][1][1]=0;
      out[25][0][1]=-210.0*tmp16*dytmp3;
      out[25][1][1]=0;
      out[26][0][1]=-350.0*tmp16*dytmp4;
      out[26][1][1]=0;
      out[27][0][1]=-490.0*tmp16*dytmp5;
      out[27][1][1]=0;
      out[28][0][1]=0;
      out[28][1][1]=10.0*dytmp17;
      out[29][0][1]=0;
      out[29][1][1]=-30.0*dytmp18;
      out[30][0][1]=0;
      out[30][1][1]=-70.0*dytmp19;
      out[31][0][1]=0;
      out[31][1][1]=-30.0*tmp9*dytmp20;
      out[32][0][1]=0;
      out[32][1][1]=-90.0*tmp9*dytmp18;
      out[33][0][1]=0;
      out[33][1][1]=-210.0*tmp9*dytmp19;
      out[34][0][1]=0;
      out[34][1][1]=-50.0*tmp10*dytmp20;
      out[35][0][1]=0;
      out[35][1][1]=-150.0*tmp10*dytmp18;
      out[36][0][1]=0;
      out[36][1][1]=-350.0*tmp10*dytmp19;
      out[37][0][1]=0;
      out[37][1][1]=-70.0*tmp11*dytmp20;
      out[38][0][1]=0;
      out[38][1][1]=-210.0*tmp11*dytmp18;
      out[39][0][1]=0;
      out[39][1][1]=-490.0*tmp11*dytmp19;

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
