// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS4_CUBE2D_LOCALBASIS_HH
#define DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS4_CUBE2D_LOCALBASIS_HH

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
  class RT4Cube2DLocalBasis
  {

  public:
    typedef LocalBasisTraits<D,2,Dune::FieldVector<D,2>,R,2,Dune::FieldVector<R,2>,
        Dune::FieldMatrix<R,2,2>, 0 > Traits;

    //! \brief Standard constructor
    RT4Cube2DLocalBasis ()
    {
      sign0 = sign1 = sign2 = sign3 = 1.0;
    }

    /**
     * \brief Make set number s, where 0 <= s < 16
     *
     * \param s Edge orientation indicator
     */
    RT4Cube2DLocalBasis (unsigned int s)
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
    constexpr unsigned int size () const
    {
      return 60;
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
      out.resize(60);

      double l1_x=(-1.0+2.0*in[0]);
      double l1_y=(-1.0+2.0*in[1]);
      double l2_x=(1.0-6.0*in[0]+6.0*pow(in[0],2));
      double l2_y=(1.0-6.0*in[1]+6.0*pow(in[1],2));
      double l3_x=(-1.0+12.0*in[0]-30.0*pow(in[0],2)+20.0*pow(in[0],3));
      double l3_y=(-1.0+12.0*in[1]-30.0*pow(in[1],2)+20.0*pow(in[1],3));
      double l4_x=(1.0-20.0*in[0]+90.0*pow(in[0],2)-140.0*pow(in[0],3)+70.0*pow(in[0],4));
      double l4_y=(1.0-20.0*in[1]+90.0*pow(in[1],2)-140.0*pow(in[1],3)+70.0*pow(in[1],4));
      double l5_x=(-1.0+30.0*in[0]-210.0*pow(in[0],2)+560.0*pow(in[0],3)-630.0*pow(in[0],4)+252.0*pow(in[0],5));
      double l5_y=(-1.0+30.0*in[1]-210.0*pow(in[1],2)+560.0*pow(in[1],3)-630.0*pow(in[1],4)+252.0*pow(in[1],5));

      out[0][0]=sign0*(0.5*(-l4_x)+0.5*l5_x);
      out[0][1]=0.0;
      out[1][0]=-(1.5)*l4_x*l1_y+1.5*l5_x*l1_y;
      out[1][1]=0.0;
      out[2][0]=sign0*(-(2.5)*l4_x*l2_y+2.5*l5_x*l2_y);
      out[2][1]=0.0;
      out[3][0]=-(3.5)*l4_x*l3_y+3.5*l5_x*l3_y;
      out[3][1]=0.0;
      out[4][0]=sign0*(-(4.5)*l4_x*l4_y+4.5*l5_x*l4_y);
      out[4][1]=0.0;

      out[5][0]=sign1*(0.5*l4_x+0.5*l5_x);
      out[5][1]=0.0;
      out[6][0]=-(1.5)*l4_x*l1_y-1.5*l5_x*l1_y;
      out[6][1]=0.0;
      out[7][0]=sign1*(2.5*l4_x*l2_y+2.5*l5_x*l2_y);
      out[7][1]=0.0;
      out[8][0]=-(3.5)*l4_x*l3_y-3.5*l5_x*l3_y;
      out[8][1]=0.0;
      out[9][0]=sign1*(4.5*l4_x*l4_y+4.5*l5_x*l4_y);
      out[9][1]=0.0;

      out[10][0]=0.0;
      out[10][1]=sign2*(0.5*(-l4_y)+0.5*l5_y);
      out[11][0]=0.0;
      out[11][1]=1.5*l1_x*l4_y-1.5*l1_x*l5_y;
      out[12][0]=0.0;
      out[12][1]=sign2*(-(2.5)*l2_x*l4_y+2.5*l2_x*l5_y);
      out[13][0]=0.0;
      out[13][1]=3.5*l3_x*l4_y-3.5*l3_x*l5_y;
      out[14][0]=0.0;
      out[14][1]=sign2*(-(4.5)*l4_x*l4_y+4.5*l4_x*l5_y);

      out[15][0]=0.0;
      out[15][1]=sign3*(0.5*l4_y+0.5*l5_y);
      out[16][0]=0.0;
      out[16][1]=1.5*l1_x*l4_y+1.5*l1_x*l5_y;
      out[17][0]=0.0;
      out[17][1]=sign3*(2.5*l2_x*l4_y+2.5*l2_x*l5_y);
      out[18][0]=0.0;
      out[18][1]=3.5*l3_x*l4_y+3.5*l3_x*l5_y;
      out[19][0]=0.0;
      out[19][1]=sign3*(4.5*l4_x*l4_y+4.5*l4_x*l5_y);

      out[20][0]=1.0-l4_x;
      out[20][1]=0.0;
      out[21][0]=3.0*l1_y-3.0*l4_x*l1_y;
      out[21][1]=0.0;
      out[22][0]=5.0*l2_y-5.0*l4_x*l2_y;
      out[22][1]=0.0;
      out[23][0]=7.0*l3_y-7.0*l4_x*l3_y;
      out[23][1]=0.0;
      out[24][0]=9.0*l4_y-9.0*l4_x*l4_y;
      out[24][1]=0.0;
      out[25][0]=3.0*l1_x-3.0*l5_x;
      out[25][1]=0.0;
      out[26][0]=9.0*l1_x*l1_y-9.0*l5_x*l1_y;
      out[26][1]=0.0;
      out[27][0]=15.0*l1_x*l2_y-15.0*l5_x*l2_y;
      out[27][1]=0.0;
      out[28][0]=21.0*l1_x*l3_y-21.0*l5_x*l3_y;
      out[28][1]=0.0;
      out[29][0]=27.0*l1_x*l4_y-27.0*l5_x*l4_y;
      out[29][1]=0.0;
      out[30][0]=5.0*l2_x-5.0*l4_x;
      out[30][1]=0.0;
      out[31][0]=15.0*l2_x*l1_y-15.0*l4_x*l1_y;
      out[31][1]=0.0;
      out[32][0]=25.0*l2_x*l2_y-25.0*l4_x*l2_y;
      out[32][1]=0.0;
      out[33][0]=35.0*l2_x*l3_y-35.0*l4_x*l3_y;
      out[33][1]=0.0;
      out[34][0]=45.0*l2_x*l4_y-45.0*l4_x*l4_y;
      out[34][1]=0.0;
      out[35][0]=7.0*l3_x-7.0*l5_x;
      out[35][1]=0.0;
      out[36][0]=21.0*l3_x*l1_y-21.0*l5_x*l1_y;
      out[36][1]=0.0;
      out[37][0]=35.0*l3_x*l2_y-35.0*l5_x*l2_y;
      out[37][1]=0.0;
      out[38][0]=49.0*l3_x*l3_y-49.0*l5_x*l3_y;
      out[38][1]=0.0;
      out[39][0]=63.0*l3_x*l4_y-63.0*l5_x*l4_y;
      out[39][1]=0.0;
      out[40][0]=0.0;
      out[40][1]=1.0-l4_y;
      out[41][0]=0.0;
      out[41][1]=3.0*l1_y-3.0*l5_y;
      out[42][0]=0.0;
      out[42][1]=5.0*l2_y-5.0*l4_y;
      out[43][0]=0.0;
      out[43][1]=7.0*l3_y-7.0*l5_y;
      out[44][0]=0.0;
      out[44][1]=3.0*l1_x-3.0*l1_x*l4_y;
      out[45][0]=0.0;
      out[45][1]=9.0*l1_x*l1_y-9.0*l1_x*l5_y;
      out[46][0]=0.0;
      out[46][1]=15.0*l1_x*l2_y-15.0*l1_x*l4_y;
      out[47][0]=0.0;
      out[47][1]=21.0*l1_x*l3_y-21.0*l1_x*l5_y;
      out[48][0]=0.0;
      out[48][1]=5.0*l2_x-5.0*l2_x*l4_y;
      out[49][0]=0.0;
      out[49][1]=15.0*l2_x*l1_y-15.0*l2_x*l5_y;
      out[50][0]=0.0;
      out[50][1]=25.0*l2_x*l2_y-25.0*l2_x*l4_y;
      out[51][0]=0.0;
      out[51][1]=35.0*l2_x*l3_y-35.0*l2_x*l5_y;
      out[52][0]=0.0;
      out[52][1]=7.0*l3_x-7.0*l3_x*l4_y;
      out[53][0]=0.0;
      out[53][1]=21.0*l3_x*l1_y-21.0*l3_x*l5_y;
      out[54][0]=0.0;
      out[54][1]=35.0*l3_x*l2_y-35.0*l3_x*l4_y;
      out[55][0]=0.0;
      out[55][1]=49.0*l3_x*l3_y-49.0*l3_x*l5_y;
      out[56][0]=0.0;
      out[56][1]=9.0*l4_x-9.0*l4_x*l4_y;
      out[57][0]=0.0;
      out[57][1]=27.0*l4_x*l1_y-27.0*l4_x*l5_y;
      out[58][0]=0.0;
      out[58][1]=45.0*l4_x*l2_y-45.0*l4_x*l4_y;
      out[59][0]=0.0;
      out[59][1]=63.0*l4_x*l3_y-63.0*l4_x*l5_y;
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
      out.resize(60);
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

    template <std::size_t dOrder>
    inline void evaluate (const std::array<int, dOrder>& /*directions*/,
                          const typename Traits::DomainType& in,         // position
                          std::vector<typename Traits::RangeType>& out) const      // return value
    {
      if (dOrder == 0) {
        evaluateFunction(in, out);
      } else {
        DUNE_THROW(NotImplemented, "Desired derivative order is not implemented");
      }
    }

    //! \brief Polynomial order of the shape functions
    unsigned int order () const
    {
      return 9;
    }

  private:
    R sign0, sign1, sign2, sign3;
  };
}

#endif // DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS3_CUBE2D_LOCALBASIS_HH
