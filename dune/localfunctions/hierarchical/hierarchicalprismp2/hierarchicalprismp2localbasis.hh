// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_HIERARCHICAL_PRISM_P2_LOCALBASIS_HH
#define DUNE_HIERARCHICAL_PRISM_P2_LOCALBASIS_HH

/** \file
    \brief Hierarchical prism p2 shape functions for the simplex
 */

#include <numeric>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <dune/localfunctions/common/localbasis.hh>

namespace Dune
{
  template<class D, class R>
  class HierarchicalPrismP2LocalBasis
  {
  public:
    //! \brief export type traits for function signature
    typedef LocalBasisTraits<D,3,Dune::FieldVector<D,3>,R,1,Dune::FieldVector<R,1>, Dune::FieldMatrix<R,1,3> > Traits;

    //! \brief number of shape functions
    unsigned int size () const
    {
      return 18;
    }

    //! \brief Evaluate all shape functions
    void evaluateFunction (const typename Traits::DomainType& in,
                           std::vector<typename Traits::RangeType> & out) const
    {
      out.resize(18);

      out[0]=(1.0-in[0]-in[1])*(1.0-in[2]);
      out[1]= in[0]*(1-in[2]);
      out[2]=in[1]*(1-in[2]);
      out[3]=in[2]*(1.0-in[0]-in[1]);
      out[4]=in[0]*in[2];
      out[5]=in[1]*in[2];

      //edges
      out[6]=2*(1.0-in[0]-in[1])*(0.5-in[0]-in[1])*(4*in[2]-4*in[2]*in[2]);
      out[7]=2*in[0]*(-0.5+in[0])*(4*in[2]-4*in[2]*in[2]);
      out[8]=2*in[1]*(-0.5+in[1])*(4*in[2]-4*in[2]*in[2]);
      out[9]=4*in[0]*(1-in[0]-in[1])*(1-3*in[2]+2*in[2]*in[2]);
      out[10]=4*in[1]*(1-in[0]-in[1])*(1-3*in[2]+2*in[2]*in[2]);
      out[11]=4*in[0]*in[1]*(1-3*in[2]+2*in[2]*in[2]);
      out[12]=4*in[0]*(1-in[0]-in[1])*(-in[2]+2*in[2]*in[2]);
      out[13]=4*in[1]*(1-in[0]-in[1])*(-in[2]+2*in[2]*in[2]);
      out[14]=4*in[0]*in[1]*(-in[2]+2*in[2]*in[2]);

      //faces
      out[15]=4*in[0]*(1-in[0]-in[1])*(4*in[2]-4*in[2]*in[2]);
      out[16]=4*in[1]*(1-in[0]-in[1])*(4*in[2]-4*in[2]*in[2]);
      out[17]=4*in[0]*in[1]*(4*in[2]-4*in[2]*in[2]);
    }



    //! \brief Evaluate Jacobian of all shape functions
    void evaluateJacobian (const typename Traits::DomainType& in,     //position
                           std::vector<typename Traits::JacobianType>& out) const  //return  value
    {
      out.resize(18);

      //vertices
      out[0][0][0] = in[2]-1;
      out[0][0][1] = in[2]-1;
      out[0][0][2] = in[0]+in[1]-1;

      out[1][0][0] = 1-in[2];
      out[1][0][1] = 0;
      out[1][0][2] =-in[0];

      out[2][0][0] = 0;
      out[2][0][1] = 1-in[2];
      out[2][0][2] = -in[1];

      out[3][0][0] = -in[2];
      out[3][0][1] = -in[2];
      out[3][0][2] = 1-in[0]-in[1];

      out[4][0][0] = in[2];
      out[4][0][1] = 0;
      out[4][0][2] = in[0];

      out[5][0][0] = 0;
      out[5][0][1] = in[2];
      out[5][0][2] = in[1];

      //edges
      out[6][0][0] = (-3+4*in[0]+4*in[1])*(4*in[2]-4*in[2]*in[2]);
      out[6][0][1] = (-3+4*in[0]+4*in[1])*(4*in[2]-4*in[2]*in[2]);
      out[6][0][2] = 2*(1-in[0]-in[1])*(0.5-in[0]-in[1])*(4-8*in[2]);

      out[7][0][0] = (-1+4*in[0])*(4*in[2]-4*in[2]*in[2]);
      out[7][0][1] = 0;
      out[7][0][2] = 2*in[0]*(-0.5+in[0])*(4-8*in[2]);

      out[8][0][0] = 0;
      out[8][0][1] = (-1+4*in[1])*(4*in[2]-4*in[2]*in[2]);
      out[8][0][2] = 2*in[1]*(-0.5+in[1])*(4-8*in[2]);

      out[9][0][0] = (4-8*in[0]-4*in[1])*(1-3*in[2]+2*in[2]*in[2]);
      out[9][0][1] = -4*in[0]*(1-3*in[2]+2*in[2]*in[2]);
      out[9][0][2] = 4*in[0]*(1-in[0]-in[1])*(-3+4*in[2]);

      out[10][0][0] = (-4*in[1])*(1-3*in[2]+2*in[2]*in[2]);
      out[10][0][1] = (4-4*in[0]-8*in[1])*(1-3*in[2]+2*in[2]*in[2]);
      out[10][0][2] = 4*in[1]*(1-in[0]-in[1])*(-3+4*in[2]);

      out[11][0][0] = 4*in[1]*(1-3*in[2]+2*in[2]*in[2]);
      out[11][0][1] = 4*in[0]*(1-3*in[2]+2*in[2]*in[2]);
      out[11][0][2] = 4*in[0]*in[1]*(-3+4*in[2]);

      out[12][0][0] = (4-8*in[0]-4*in[1])*(-in[2]+2*in[2]*in[2]);
      out[12][0][1] = (-4*in[0])*(-in[2]+2*in[2]*in[2]);
      out[12][0][2] = 4*in[0]*(1-in[0]-in[1])*(-1+4*in[2]);

      out[13][0][0] = -4*in[1]*(-in[2]+2*in[2]*in[2]);
      out[13][0][1] = (4-4*in[0]-8*in[1])*(-in[2]+2*in[2]*in[2]);
      out[13][0][2] = 4*in[1]*(1-in[0]-in[1])*(-1+4*in[2]);

      out[14][0][0] = 4*in[1]*(-in[2]+2*in[2]*in[2]);
      out[14][0][1] = 4*in[0]*(-in[2]+2*in[2]*in[2]);
      out[14][0][2] = 4*in[0]*in[1]*(-1+4*in[2]);

      //faces
      out[15][0][0] = (4-8*in[0]-4*in[1])*(4*in[2]-4*in[2]*in[2]);
      out[15][0][1] = -4*in[0]*(4*in[2]-4*in[2]*in[2]);
      out[15][0][2] = 4*in[0]*(1-in[0]-in[1])*(4-8*in[2]);

      out[16][0][0] = -4*in[1]*(4*in[2]-4*in[2]*in[2]);
      out[16][0][1] = (4-4*in[0]-8*in[1])*(4*in[2]-4*in[2]*in[2]);
      out[16][0][2] = 4*in[1]*(1-in[0]-in[1])*(4-8*in[2]);

      out[17][0][0] = 4*in[1]*(4*in[2]-4*in[2]*in[2]);
      out[17][0][1] = 4*in[0]*(4*in[2]-4*in[2]*in[2]);
      out[17][0][2] = 4*in[0]*in[1]*(4-8*in[2]);
    }

    //! \brief Evaluate partial derivatives of all shape functions
    void partial (const std::array<unsigned int, 3>& order,
                  const typename Traits::DomainType& in,         // position
                  std::vector<typename Traits::RangeType>& out) const      // return value
    {
      auto totalOrder = std::accumulate(order.begin(), order.end(), 0);
      if (totalOrder == 0) {
        evaluateFunction(in, out);
      } else if (totalOrder == 1) {
        out.resize(size());
        auto const direction = std::distance(order.begin(), std::find(order.begin(), order.end(), 1));

        switch (direction) {
        case 0:
          out[0]  = in[2]-1;
          out[1]  = 1-in[2];
          out[2]  = 0;
          out[3]  = -in[2];
          out[4]  = in[2];
          out[5]  = 0;
          out[6]  = (-3+4*in[0]+4*in[1])*(4*in[2]-4*in[2]*in[2]);
          out[7]  = (-1+4*in[0])*(4*in[2]-4*in[2]*in[2]);
          out[8]  = 0;
          out[9]  = (4-8*in[0]-4*in[1])*(1-3*in[2]+2*in[2]*in[2]);
          out[10] = (-4*in[1])*(1-3*in[2]+2*in[2]*in[2]);
          out[11] = 4*in[1]*(1-3*in[2]+2*in[2]*in[2]);
          out[12] = (4-8*in[0]-4*in[1])*(-in[2]+2*in[2]*in[2]);
          out[13] = -4*in[1]*(-in[2]+2*in[2]*in[2]);
          out[14] = 4*in[1]*(-in[2]+2*in[2]*in[2]);
          out[15] = (4-8*in[0]-4*in[1])*(4*in[2]-4*in[2]*in[2]);
          out[16] = -4*in[1]*(4*in[2]-4*in[2]*in[2]);
          out[17] = 4*in[1]*(4*in[2]-4*in[2]*in[2]);
          break;
        case 1:
          out[0]  = in[2]-1;
          out[1]  = 0;
          out[2]  = 1-in[2];
          out[3]  = -in[2];
          out[4]  = 0;
          out[5]  = in[2];
          out[6]  = (-3+4*in[0]+4*in[1])*(4*in[2]-4*in[2]*in[2]);
          out[7]  = 0;
          out[8]  = (-1+4*in[1])*(4*in[2]-4*in[2]*in[2]);
          out[9]  = -4*in[0]*(1-3*in[2]+2*in[2]*in[2]);
          out[10] = (4-4*in[0]-8*in[1])*(1-3*in[2]+2*in[2]*in[2]);
          out[11] = 4*in[0]*(1-3*in[2]+2*in[2]*in[2]);
          out[12] = (-4*in[0])*(-in[2]+2*in[2]*in[2]);
          out[13] = (4-4*in[0]-8*in[1])*(-in[2]+2*in[2]*in[2]);
          out[14] = 4*in[0]*(-in[2]+2*in[2]*in[2]);
          out[15] = -4*in[0]*(4*in[2]-4*in[2]*in[2]);
          out[16] = (4-4*in[0]-8*in[1])*(4*in[2]-4*in[2]*in[2]);
          out[17] = 4*in[0]*(4*in[2]-4*in[2]*in[2]);
          break;
        case 2:
          out[0]  = in[0]+in[1]-1;
          out[1]  =-in[0];
          out[2]  = -in[1];
          out[3]  = 1-in[0]-in[1];
          out[4]  = in[0];
          out[5]  = in[1];
          out[6]  = 2*(1-in[0]-in[1])*(0.5-in[0]-in[1])*(4-8*in[2]);
          out[7]  = 2*in[0]*(-0.5+in[0])*(4-8*in[2]);
          out[8]  = 2*in[1]*(-0.5+in[1])*(4-8*in[2]);
          out[9]  = 4*in[0]*(1-in[0]-in[1])*(-3+4*in[2]);
          out[10] = 4*in[1]*(1-in[0]-in[1])*(-3+4*in[2]);
          out[11] = 4*in[0]*in[1]*(-3+4*in[2]);
          out[12] = 4*in[0]*(1-in[0]-in[1])*(-1+4*in[2]);
          out[13] = 4*in[1]*(1-in[0]-in[1])*(-1+4*in[2]);
          out[14] = 4*in[0]*in[1]*(-1+4*in[2]);
          out[15] = 4*in[0]*(1-in[0]-in[1])*(4-8*in[2]);
          out[16] = 4*in[1]*(1-in[0]-in[1])*(4-8*in[2]);
          out[17] = 4*in[0]*in[1]*(4-8*in[2]);
          break;
        default:
          DUNE_THROW(RangeError, "Component out of range.");
        }
      } else {
        DUNE_THROW(NotImplemented, "Desired derivative order is not implemented");
      }
    }

    /** \brief Polynomial order of the shape functions
     */
    unsigned int order() const
    {
      return 2;
    }

  };
}
#endif
