// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS0_CUBE3D_ALL_HH
#define DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS0_CUBE3D_ALL_HH

#include <cstddef>
#include <numeric>
#include <vector>

#include <dune/common/fmatrix.hh>

#include <dune/localfunctions/common/localbasis.hh>
#include <dune/localfunctions/common/localkey.hh>
#include <dune/localfunctions/common/localinterpolation.hh>

namespace Dune
{
  /**@ingroup LocalBasisImplementation
         \brief Lowest order Raviart-Thomas shape functions on the reference hexahedron.

         \tparam D Type to represent the field in the domain.
         \tparam R Type to represent the field in the range.

         \nosubgrouping
   */
  template<class D, class R>
  class RT0Cube3DLocalBasis
  {
  public:
    typedef LocalBasisTraits<D,3,Dune::FieldVector<D,3>,R,3,Dune::FieldVector<R,3>,
        Dune::FieldMatrix<R,3,3> > Traits;

    //! \brief Make set number s, where 0 <= s < 64
    RT0Cube3DLocalBasis (unsigned int s = 0)
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

    //! \brief Evaluate partial derivatives of all shape functions
    void partial (const std::array<unsigned int, 3>& order,
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
          out[i][0] = out[i][1] = out[i][2] = 0;

        switch (direction) {
        case 0:
          out[0][0] = sign0;
          out[1][0] = sign1;
          break;
        case 1:
          out[2][1] = sign2;
          out[3][1] = sign3;
          break;
        case 2:
          out[4][2] = sign4;
          out[5][2] = sign5;
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
    R sign0, sign1, sign2, sign3, sign4, sign5;
  };


  /**@ingroup LocalInterpolationImplementation
         \brief Lowest order Raviart-Thomas shape functions on the reference hexahedron.

         \tparam LB corresponding LocalBasis giving traits

         \nosubgrouping
   */
  template<class LB>
  class RT0Cube3DLocalInterpolation
  {
  public:

    //! \brief Make set number s, where 0 <= s < 64
    RT0Cube3DLocalInterpolation (unsigned int s = 0)
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
    void interpolate (const F& ff, std::vector<C>& out) const
    {
      // f gives v*outer normal at a point on the edge!
      auto&& f = Impl::makeFunctionWithCallOperator<typename LB::Traits::DomainType>(ff);

      out.resize(6);

      auto y = f(m0); out[0] = (y[0]*n0[0]+y[1]*n0[1]+y[2]*n0[2])*sign0;
      y = f(m1); out[1] = (y[0]*n1[0]+y[1]*n1[1]+y[2]*n1[2])*sign1;
      y = f(m2); out[2] = (y[0]*n2[0]+y[1]*n2[1]+y[2]*n2[2])*sign2;
      y = f(m3); out[3] = (y[0]*n3[0]+y[1]*n3[1]+y[2]*n3[2])*sign3;
      y = f(m4); out[4] = (y[0]*n4[0]+y[1]*n4[1]+y[2]*n4[2])*sign4;
      y = f(m5); out[5] = (y[0]*n5[0]+y[1]*n5[1]+y[2]*n5[2])*sign5;
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
  class RT0Cube3DLocalCoefficients
  {
  public:
    //! \brief Standard constructor
    RT0Cube3DLocalCoefficients () : li(6)
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
#endif // DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS0_CUBE3D_ALL_HH
