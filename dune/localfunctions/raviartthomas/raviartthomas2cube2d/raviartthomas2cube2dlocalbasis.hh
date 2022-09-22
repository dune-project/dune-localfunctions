// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS2_CUBE2D_LOCALBASIS_HH
#define DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS2_CUBE2D_LOCALBASIS_HH

#include <numeric>
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
  class RT2Cube2DLocalBasis
  {

  public:
    typedef LocalBasisTraits<D,2,Dune::FieldVector<D,2>,R,2,Dune::FieldVector<R,2>,
        Dune::FieldMatrix<R,2,2> > Traits;

    /**
     * \brief Make set number s, where 0 <= s < 16
     *
     * \param s Edge orientation indicator
     */
    RT2Cube2DLocalBasis (unsigned int s = 0)
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
      return 24;
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
      out.resize(24);

      out[0][0] = sign0*(-1.0 + 9.0*in[0] - 18.0*in[0]*in[0] + 10.0*in[0]*in[0]*in[0]);
      out[0][1] = 0.0;
      out[1][0] = 3.0 - 27.0*in[0] - 6.0*in[1] + 54.0*in[0]*in[1] + 54.0*in[0]*in[0] - 108.0*in[0]*in[0]*in[1] - 30.0*in[0]*in[0]*in[0] + 60.0*in[0]*in[0]*in[0]*in[1];
      out[1][1] = 0.0;
      out[2][0] = sign0*(-5.0 + 45.0*in[0] + 30.0*in[1] - 270.0*in[0]*in[1] - 90.0*in[0]*in[0] - 30.0*in[1]*in[1] + 540.0*in[0]*in[0]*in[1] + 270.0*in[0]*in[1]*in[1] + 50.0*in[0]*in[0]*in[0] - 540.0*in[0]*in[0]*in[1]*in[1] - 300.0*in[0]*in[0]*in[0]*in[1] + 300.0*in[0]*in[0]*in[0]*in[1]*in[1]);
      out[2][1] = 0.0;
      out[3][0] = sign1*(3.0*in[0] - 12.0*in[0]*in[0] + 10.0*in[0]*in[0]*in[0]);
      out[3][1] = 0.0;
      out[4][0] = 9.0*in[0] - 18.0*in[0]*in[1] - 36.0*in[0]*in[0] + 72.0*in[0]*in[0]*in[1] + 30.0*in[0]*in[0]*in[0] - 60.0*in[0]*in[0]*in[0]*in[1];
      out[4][1] = 0.0;
      out[5][0] = sign1*(15.0*in[0] - 90.0*in[0]*in[1] - 60.0*in[0]*in[0] + 360.0*in[0]*in[0]*in[1] + 90.0*in[0]*in[1]*in[1] + 50.0*in[0]*in[0]*in[0] - 360.0*in[0]*in[0]*in[1]*in[1] - 300.0*in[0]*in[0]*in[0]*in[1] + 300.0*in[0]*in[0]*in[0]*in[1]*in[1]);
      out[5][1] = 0.0;
      out[6][0] = 0.0;
      out[6][1] = sign2*(-1.0 + 9.0*in[1] - 18.0*in[1]*in[1] + 10.0*in[1]*in[1]*in[1]);
      out[7][0] = 0.0;
      out[7][1] = -3.0 + 6.0*in[0] + 27.0*in[1] - 54.0*in[0]*in[1] - 54.0*in[1]*in[1] + 108.0*in[0]*in[1]*in[1] + 30.0*in[1]*in[1]*in[1] - 60.0*in[0]*in[1]*in[1]*in[1];
      out[8][0] = 0.0;
      out[8][1] = sign2*(-5.0 + 30.0*in[0] + 45.0*in[1] - 270.0*in[0]*in[1] - 30.0*in[0]*in[0] - 90.0*in[1]*in[1] + 270.0*in[0]*in[0]*in[1] + 540.0*in[0]*in[1]*in[1] + 50.0*in[1]*in[1]*in[1] - 540.0*in[0]*in[0]*in[1]*in[1] - 300.0*in[0]*in[1]*in[1]*in[1] + 300.0*in[1]*in[1]*in[1]*in[0]*in[0]);
      out[9][0] = 0.0;
      out[9][1] = sign3*(3.0*in[1] - 12.0*in[1]*in[1] + 10.0*in[1]*in[1]*in[1]);
      out[10][0] = 0.0;
      out[10][1] = -9.0*in[1] + 18.0*in[0]*in[1] + 36.0*in[1]*in[1] - 72.0*in[0]*in[1]*in[1] - 30.0*in[1]*in[1]*in[1] + 60.0*in[0]*in[1]*in[1]*in[1];
      out[11][0] = 0.0;
      out[11][1] = sign3*(15.0*in[1] - 90.0*in[0]*in[1] - 60.0*in[1]*in[1] + 90.0*in[0]*in[0]*in[1] + 360.0*in[0]*in[1]*in[1] + 50.0*in[1]*in[1]*in[1] - 360.0*in[0]*in[0]*in[1]*in[1] - 300.0*in[0]*in[1]*in[1]*in[1] + 300.0*in[1]*in[1]*in[1]*in[0]*in[0]);
      out[12][0] = 324.0*in[0] -1296.0*in[0]*in[1] - 864.0*in[0]*in[0] + 3456.0*in[0]*in[0]*in[1] + 1080.0*in[0]*in[1]*in[1] + 540.0*in[0]*in[0]*in[0] - 2880.0*in[0]*in[0]*in[1]*in[1] - 2160.0*in[0]*in[0]*in[0]*in[1] + 1800.0*in[0]*in[0]*in[0]*in[1]*in[1];
      out[12][1] = 0.0;
      out[13][0] = 0.0;
      out[13][1] = 324.0*in[1] - 1296.0*in[0]*in[1] - 864.0*in[1]*in[1] + 1080.0*in[0]*in[0]*in[1] + 3456.0*in[0]*in[1]*in[1] + 540.0*in[1]*in[1]*in[1] - 2880.0*in[0]*in[0]*in[1]*in[1] - 2160.0*in[0]*in[1]*in[1]*in[1] + 1800.0*in[1]*in[1]*in[1]*in[0]*in[0];
      out[14][0] = -540.0*in[0] + 2160.0*in[0]*in[1] + 1620.0*in[0]*in[0] - 6480.0*in[0]*in[0]*in[1] - 1800.0*in[0]*in[1]*in[1] - 1080.0*in[0]*in[0]*in[0] + 5400.0*in[0]*in[0]*in[1]*in[1] + 4320.0*in[0]*in[0]*in[0]*in[1] - 3600.0*in[0]*in[0]*in[0]*in[1]*in[1];
      out[14][1] = 0.0;
      out[15][0] = 0.0;
      out[15][1] = -1296.0*in[1] + 6912.0*in[0]*in[1] + 3456.0*in[1]*in[1] - 6480.0*in[0]*in[0]*in[1] - 18432.0*in[0]*in[1]*in[1] - 2160.0*in[1]*in[1]*in[1] + 17280.0*in[0]*in[0]*in[1]*in[1] + 11520.0*in[0]*in[1]*in[1]*in[1] - 10800.0*in[0]*in[0]*in[1]*in[1]*in[1];
      out[16][0] = -1296.0*in[0] + 6912.0*in[0]*in[1] + 3456.0*in[0]*in[0] - 6480.0*in[0]*in[1]*in[1] - 18432.0*in[0]*in[0]*in[1] - 2160.0*in[0]*in[0]*in[0] + 17280.0*in[0]*in[0]*in[1]*in[1] + 11520.0*in[1]*in[0]*in[0]*in[0] - 10800.0*in[0]*in[0]*in[0]*in[1]*in[1];
      out[16][1] = 0.0;
      out[17][0] = 0.0;
      out[17][1] = -540.0*in[1] + 2160.0*in[0]*in[1] + 1620.0*in[1]*in[1] - 1800.0*in[0]*in[0]*in[1] - 6480.0*in[0]*in[1]*in[1] - 1080.0*in[1]*in[1]*in[1] + 5400.0*in[0]*in[0]*in[1]*in[1] + 4320.0*in[0]*in[1]*in[1]*in[1] - 3600.0*in[0]*in[0]*in[1]*in[1]*in[1];
      out[18][0] = 2160.0*in[0] - 11520.0*in[0]*in[1] - 6480.0*in[0]*in[0] + 34560.0*in[0]*in[0]*in[1] + 10800.0*in[0]*in[1]*in[1] + 4320.0*in[0]*in[0]*in[0] - 32400.0*in[0]*in[0]*in[1]*in[1] - 23040.0*in[0]*in[0]*in[0]*in[1] + 21600.0*in[0]*in[0]*in[0]*in[1]*in[1];
      out[18][1] = 0.0;
      out[19][0] = 0.0;
      out[19][1] = 2160.0*in[1] - 11520.0*in[0]*in[1] - 6480.0*in[1]*in[1] + 10800.0*in[0]*in[0]*in[1] + 34560.0*in[0]*in[1]*in[1] + 4320.0*in[1]*in[1]*in[1] - 32400.0*in[0]*in[0]*in[1]*in[1] - 23040.0*in[0]*in[1]*in[1]*in[1] + 21600.0*in[0]*in[0]*in[1]*in[1]*in[1];
      out[20][0] = 1080.0*in[0] - 6480.0*in[0]*in[1] - 2880.0*in[0]*in[0] + 17280.0*in[0]*in[0]*in[1] + 6480.0*in[0]*in[1]*in[1] + 1800.0*in[0]*in[0]*in[0] - 17280.0*in[0]*in[0]*in[1]*in[1] - 10800.0*in[0]*in[0]*in[0]*in[1] + 10800.0*in[0]*in[0]*in[0]*in[1]*in[1];
      out[20][1] = 0.0;
      out[21][0] = 0.0;
      out[21][1] = 1080.0*in[1] - 6480.0*in[0]*in[1] - 2880.0*in[1]*in[1] + 6480.0*in[0]*in[0]*in[1] + 17280.0*in[0]*in[1]*in[1] + 1800.0*in[1]*in[1]*in[1] - 17280.0*in[0]*in[0]*in[1]*in[1] - 10800.0*in[0]*in[1]*in[1]*in[1] + 10800.0*in[0]*in[0]*in[1]*in[1]*in[1];
      out[22][0] = -1800.0*in[0] + 10800.0*in[0]*in[1] + 5400.0*in[0]*in[0] - 32400.0*in[0]*in[0]*in[1] - 10800.0*in[0]*in[1]*in[1] - 3600.0*in[0]*in[0]*in[0] + 32400.0*in[0]*in[0]*in[1]*in[1] + 21600.0*in[0]*in[0]*in[0]*in[1] - 21600.0*in[0]*in[0]*in[0]*in[1]*in[1];
      out[22][1] = 0.0;
      out[23][0] = 0.0;
      out[23][1] = -1800.0*in[1] + 10800.0*in[0]*in[1] + 5400.0*in[1]*in[1] - 10800.0*in[0]*in[0]*in[1] - 32400.0*in[0]*in[1]*in[1] - 3600.0*in[1]*in[1]*in[1] + 32400.0*in[0]*in[0]*in[1]*in[1] + 21600.0*in[0]*in[1]*in[1]*in[1] - 21600.0*in[0]*in[0]*in[1]*in[1]*in[1];
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
      out.resize(24);

      out[0][0][0] = sign0*(9.0 - 36.0*in[0] + 30.0*in[0]*in[0]);
      out[0][0][1] = 0.0;
      out[0][1][0] = 0.0;
      out[0][1][1] = 0.0;

      out[1][0][0] = -27.0 + 54.0*in[1] + 108.0*in[0] - 216.0*in[0]*in[1] - 90.0*in[0]*in[0] + 180.0*in[0]*in[0]*in[1];
      out[1][0][1] = -6.0 + 54.0*in[0] - 108.0*in[0]*in[0] + 60.0*in[0]*in[0]*in[0];
      out[1][1][0] = 0.0;
      out[1][1][1] = 0.0;

      out[2][0][0] = sign0*(45.0 - 270.0*in[1] - 180.0*in[0] + 1080.0*in[0]*in[1] + 270.0*in[1]*in[1] + 150.0*in[0]*in[0] - 1080.0*in[0]*in[1]*in[1] - 900.0*in[0]*in[0]*in[1] + 900.0*in[0]*in[0]*in[1]*in[1]);
      out[2][0][1] = sign0*(30.0 - 270.0*in[0] - 60.0*in[1] + 540.0*in[0]*in[0] + 540.0*in[0]*in[1] - 1080.0*in[0]*in[0]*in[1] - 300.0*in[0]*in[0]*in[0] + 600.0*in[0]*in[0]*in[0]*in[1]);
      out[2][1][0] = 0.0;
      out[2][1][1] = 0.0;

      out[3][0][0] = sign1*(3.0 - 24.0*in[0] + 30.0*in[0]*in[0]);
      out[3][0][1] = 0.0;
      out[3][1][0] = 0.0;
      out[3][1][1] = 0.0;

      out[4][0][0] = 9.0 - 18.0*in[1] - 72.0*in[0] + 144.0*in[0]*in[1] + 90.0*in[0]*in[0] - 180.0*in[0]*in[0]*in[1];
      out[4][0][1] = -18.0*in[0] + 72.0*in[0]*in[0] - 60.0*in[0]*in[0]*in[0];
      out[4][1][0] = 0.0;
      out[4][1][1] = 0.0;

      out[5][0][0] = sign1*(15.0 - 90.0*in[1] - 120.0*in[0] + 720.0*in[0]*in[1] + 90.0*in[1]*in[1] + 150.0*in[0]*in[0] - 720.0*in[0]*in[1]*in[1] - 900.0*in[0]*in[0]*in[1] + 900.0*in[0]*in[0]*in[1]*in[1]);
      out[5][0][1] = sign1*(-90.0*in[0] + 360.0*in[0]*in[0] + 180.0*in[0]*in[1] - 720.0*in[0]*in[0]*in[1] - 300.0*in[0]*in[0]*in[0] + 600.0*in[0]*in[0]*in[0]*in[1]);
      out[5][1][0] = 0.0;
      out[5][1][1] = 0.0;


      out[6][0][0] = 0.0;
      out[6][0][1] = 0.0;
      out[6][1][0] = 0.0;
      out[6][1][1] = sign2*(9.0 - 36.0*in[1] + 30.0*in[1]*in[1]);

      out[7][0][0] = 0.0;
      out[7][0][1] = 0.0;
      out[7][1][0] = 6.0 - 54.0*in[1] + 108.0*in[1]*in[1] - 60.0*in[1]*in[1]*in[1];
      out[7][1][1] = 27.0 - 54.0*in[0] - 108.0*in[1] + 216.0*in[0]*in[1] + 90.0*in[1]*in[1] - 180.0*in[0]*in[1]*in[1];

      out[8][0][0] = 0.0;
      out[8][0][1] = 0.0;
      out[8][1][0] = sign2*(30.0 - 270.0*in[1] - 60.0*in[0] + 540.0*in[0]*in[1] + 540.0*in[1]*in[1] - 1080.0*in[0]*in[1]*in[1] - 300.0*in[1]*in[1]*in[1] + 600.0*in[1]*in[1]*in[1]*in[0]);
      out[8][1][1] = sign2*(45.0 - 270.0*in[0] - 180.0*in[1] + 270.0*in[0]*in[0] + 1080.0*in[0]*in[1] + 150.0*in[1]*in[1] - 1080.0*in[0]*in[0]*in[1] - 900.0*in[0]*in[1]*in[1] + 900.0*in[1]*in[1]*in[0]*in[0]);

      out[9][0][0] = 0.0;
      out[9][0][1] = 0.0;
      out[9][1][0] = 0.0;
      out[9][1][1] = sign3*(3.0 - 24.0*in[1] + 30.0*in[1]*in[1]);

      out[10][0][0] = 0.0;
      out[10][0][1] = 0.0;
      out[10][1][0] = 18.0*in[1] - 72.0*in[1]*in[1] + 60.0*in[1]*in[1]*in[1];
      out[10][1][1] = -9.0 + 18.0*in[0] + 72.0*in[1] - 144.0*in[0]*in[1] - 90.0*in[1]*in[1] + 180.0*in[0]*in[1]*in[1];

      out[11][0][0] = 0.0;
      out[11][0][1] = 0.0;
      out[11][1][0] = sign3*(-90.0*in[1] + 180.0*in[0]*in[1] + 360.0*in[1]*in[1] - 720.0*in[0]*in[1]*in[1] - 300.0*in[1]*in[1]*in[1] + 600.0*in[1]*in[1]*in[1]*in[0]);
      out[11][1][1] = sign3*(15.0 - 90.0*in[0] - 120.0*in[1] + 90.0*in[0]*in[0] + 720.0*in[0]*in[1] + 150.0*in[1]*in[1] - 720.0*in[0]*in[0]*in[1] - 900.0*in[0]*in[1]*in[1] + 900.0*in[1]*in[1]*in[0]*in[0]);

      out[12][0][0] = 324 -1296.0*in[1] - 1728.0*in[0] + 6912.0*in[0]*in[1] + 1080.0*in[1]*in[1] + 1620.0*in[0]*in[0] - 5760.0*in[0]*in[1]*in[1] - 6480.0*in[0]*in[0]*in[1] + 5400.0*in[0]*in[0]*in[1]*in[1];
      out[12][0][1] = -1296.0*in[0] + 3456.0*in[0]*in[0] + 2160.0*in[0]*in[1] - 5760.0*in[0]*in[0]*in[1] - 2160.0*in[0]*in[0]*in[0] + 3600.0*in[0]*in[0]*in[0]*in[1];
      out[12][1][0] = 0.0;
      out[12][1][1] = 0.0;

      out[13][0][0] = 0.0;
      out[13][0][1] = 0.0;
      out[13][1][0] = -1296.0*in[1] + 2160.0*in[0]*in[1] + 3456.0*in[1]*in[1] - 5760.0*in[0]*in[1]*in[1] - 2160.0*in[1]*in[1]*in[1] + 3600.0*in[1]*in[1]*in[1]*in[0];
      out[13][1][1] = 324.0 - 1296.0*in[0] - 1728.0*in[1] + 1080.0*in[0]*in[0] + 6912.0*in[0]*in[1] + 1620.0*in[1]*in[1] - 5760.0*in[0]*in[0]*in[1] - 6480.0*in[0]*in[1]*in[1] + 5400.0*in[1]*in[1]*in[0]*in[0];

      out[14][0][0] = -540.0 + 2160.0*in[1] + 3240.0*in[0] - 12960.0*in[0]*in[1] - 1800.0*in[1]*in[1] - 3240.0*in[0]*in[0] + 10800.0*in[0]*in[1]*in[1] + 12960.0*in[0]*in[0]*in[1] - 10800.0*in[0]*in[0]*in[1]*in[1];
      out[14][0][1] = 2160.0*in[0] - 6480.0*in[0]*in[0] - 3600.0*in[0]*in[1] + 10800.0*in[0]*in[0]*in[1] + 4320.0*in[0]*in[0]*in[0] - 7200.0*in[0]*in[0]*in[0]*in[1];
      out[14][1][0] = 0.0;
      out[14][1][1] = 0.0;

      out[15][0][0] = 0.0;
      out[15][0][1] = 0.0;
      out[15][1][0] = 6912.0*in[1] - 12960.0*in[0]*in[1] - 18432.0*in[1]*in[1] + 34560.0*in[0]*in[1]*in[1] + 11520.0*in[1]*in[1]*in[1] - 21600.0*in[0]*in[1]*in[1]*in[1];
      out[15][1][1] = -1296.0 + 6912.0*in[0] + 6912.0*in[1] - 6480.0*in[0]*in[0] - 36864.0*in[0]*in[1] - 6480.0*in[1]*in[1] + 34560.0*in[0]*in[0]*in[1] + 34560.0*in[0]*in[1]*in[1] - 32400.0*in[0]*in[0]*in[1]*in[1];

      out[16][0][0] = -1296.0 + 6912.0*in[1] + 6912.0*in[0] - 6480.0*in[1]*in[1] - 36864.0*in[0]*in[1] - 6480.0*in[0]*in[0] + 34560.0*in[0]*in[1]*in[1] + 34560.0*in[1]*in[0]*in[0] - 32400.0*in[0]*in[0]*in[1]*in[1];
      out[16][0][1] = 6912.0*in[0] - 12960.0*in[0]*in[1] - 18432.0*in[0]*in[0] + 34560.0*in[0]*in[0]*in[1] + 11520.0*in[0]*in[0]*in[0] - 21600.0*in[0]*in[0]*in[0]*in[1];
      out[16][1][0] = 0.0;
      out[16][1][1] = 0.0;

      out[17][0][0] = 0.0;
      out[17][0][1] = 0.0;
      out[17][1][0] = 2160.0*in[1] - 3600.0*in[0]*in[1] - 6480.0*in[1]*in[1] + 10800.0*in[0]*in[1]*in[1] + 4320.0*in[1]*in[1]*in[1] - 7200.0*in[0]*in[1]*in[1]*in[1];
      out[17][1][1] = -540.0 + 2160.0*in[0] + 3240.0*in[1] - 1800.0*in[0]*in[0] - 12960.0*in[0]*in[1] - 3240.0*in[1]*in[1] + 10800.0*in[0]*in[0]*in[1] + 12960.0*in[0]*in[1]*in[1] - 10800.0*in[0]*in[0]*in[1]*in[1];

      out[18][0][0] = 2160.0 - 11520.0*in[1] - 12960.0*in[0] + 69120.0*in[0]*in[1] + 10800.0*in[1]*in[1] + 12960.0*in[0]*in[0] - 64800.0*in[0]*in[1]*in[1] - 69120.0*in[0]*in[0]*in[1] + 64800.0*in[0]*in[0]*in[1]*in[1];
      out[18][0][1] = -11520.0*in[0] + 34560.0*in[0]*in[0] + 21600.0*in[0]*in[1] - 64800.0*in[0]*in[0]*in[1] - 23040.0*in[0]*in[0]*in[0] + 43200.0*in[0]*in[0]*in[0]*in[1];
      out[18][1][0] = 0.0;
      out[18][1][1] = 0.0;

      out[19][0][0] = 0.0;
      out[19][0][1] = 0.0;
      out[19][1][0] = -11520.0*in[1] + 21600.0*in[0]*in[1] + 34560.0*in[1]*in[1] - 64800.0*in[0]*in[1]*in[1] - 23040.0*in[1]*in[1]*in[1] + 43200.0*in[0]*in[1]*in[1]*in[1];
      out[19][1][1] = 2160.0 - 11520.0*in[0] - 12960.0*in[1] + 10800.0*in[0]*in[0] + 69120.0*in[0]*in[1] + 12960.0*in[1]*in[1] - 64800.0*in[0]*in[0]*in[1] - 69120.0*in[0]*in[1]*in[1] + 64800.0*in[0]*in[0]*in[1]*in[1];

      out[20][0][0] = 1080.0 - 6480.0*in[1] - 5760.0*in[0] + 34560.0*in[0]*in[1] + 6480.0*in[1]*in[1] + 5400.0*in[0]*in[0] - 34560.0*in[0]*in[1]*in[1] - 32400.0*in[0]*in[0]*in[1] + 32400.0*in[0]*in[0]*in[1]*in[1];
      out[20][0][1] = -6480.0*in[0] + 17280.0*in[0]*in[0] + 12960.0*in[0]*in[1] - 34560.0*in[0]*in[0]*in[1] - 10800.0*in[0]*in[0]*in[0] + 21600.0*in[0]*in[0]*in[0]*in[1];
      out[20][1][0] = 0.0;
      out[20][1][1] = 0.0;

      out[21][0][0] = 0.0;
      out[21][0][1] = 0.0;
      out[21][1][0] = -6480.0*in[1] + 12960.0*in[0]*in[1] + 17280.0*in[1]*in[1] - 34560.0*in[0]*in[1]*in[1] - 10800.0*in[1]*in[1]*in[1] + 21600.0*in[0]*in[1]*in[1]*in[1];
      out[21][1][1] = 1080.0 - 6480.0*in[0] - 5760.0*in[1] + 6480.0*in[0]*in[0] + 34560.0*in[0]*in[1] + 5400.0*in[1]*in[1] - 34560.0*in[0]*in[0]*in[1] - 32400.0*in[0]*in[1]*in[1] + 32400.0*in[0]*in[0]*in[1]*in[1];

      out[22][0][0] = -1800.0 + 10800.0*in[1] + 10800.0*in[0] - 64800.0*in[0]*in[1] - 10800.0*in[1]*in[1] - 10800.0*in[0]*in[0] + 64800.0*in[0]*in[1]*in[1] + 64800.0*in[0]*in[0]*in[1] - 64800.0*in[0]*in[0]*in[1]*in[1];
      out[22][0][1] = 10800.0*in[0] - 32400.0*in[0]*in[0] - 21600.0*in[0]*in[1] + 64800.0*in[0]*in[0]*in[1] + 21600.0*in[0]*in[0]*in[0] - 43200.0*in[0]*in[0]*in[0]*in[1];
      out[22][1][0] = 0.0;
      out[22][1][1] = 0.0;

      out[23][0][0] = 0.0;
      out[23][0][1] = 0.0;
      out[23][1][0] = 10800.0*in[1] - 21600.0*in[0]*in[1] - 32400.0*in[1]*in[1] + 64800.0*in[0]*in[1]*in[1] + 21600.0*in[1]*in[1]*in[1] - 43200.0*in[0]*in[1]*in[1]*in[1];
      out[23][1][1] = -1800.0 + 10800.0*in[0] + 10800.0*in[1] - 10800.0*in[0]*in[0] - 64800.0*in[0]*in[1] - 10800.0*in[1]*in[1] + 64800.0*in[0]*in[0]*in[1] + 64800.0*in[0]*in[1]*in[1] - 64800.0*in[0]*in[0]*in[1]*in[1];
    }

    //! \brief Evaluate partial derivatives of all shape functions
    void partial (const std::array<unsigned int, 2>& order,
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

    //! \brief Polynomial order of the shape functions
    unsigned int order () const
    {
      return 5;
    }

  private:
    R sign0, sign1, sign2, sign3;
  };
}
#endif // DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS2_CUBE2D_LOCALBASIS_HH
