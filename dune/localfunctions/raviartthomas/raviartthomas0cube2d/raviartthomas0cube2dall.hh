// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS0_CUBE2D_ALL_HH
#define DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS0_CUBE2D_ALL_HH

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

    //! \brief Constructor with a set of edge orientations
    RT0Cube2DLocalBasis (std::bitset<4> s = 0)
    {
      for (int i=0; i<4; i++)
        sign_[i] = s[i] ? -1.0 : 1.0;
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
      out[0] = {sign_[0]*(in[0]-1.0), 0.0};
      out[1] = {sign_[1]*(in[0]),     0.0};
      out[2] = {0.0,                  sign_[2]*(in[1]-1.0)};
      out[3] = {0.0,                  sign_[3]*(in[1])};
    }

    //! \brief Evaluate Jacobian of all shape functions
    inline void
    evaluateJacobian (const typename Traits::DomainType& in,             // position
                      std::vector<typename Traits::JacobianType>& out) const                          // return value
    {
      out.resize(4);
      out[0][0] = {sign_[0], 0};
      out[0][1] = {0,        0};

      out[1][0] = {sign_[1], 0};
      out[1][1] = {0,        0};

      out[2][0] = {0,        0};
      out[2][1] = {0, sign_[2]};

      out[3][0] = {0,        0};
      out[3][1] = {0, sign_[3]};
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
          out[i] = {0, 0};

        switch (direction) {
        case 0:
          out[0][0] = sign_[0];
          out[1][0] = sign_[1];
          break;
        case 1:
          out[2][1] = sign_[2];
          out[3][1] = sign_[3];
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
    std::array<R,4> sign_;
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

    //! \brief Constructor with explicitly given edge orientations
    RT0Cube2DLocalInterpolation (std::bitset<4> s = 0)
    {
      for (int i=0; i<4; i++)
        sign_[i] = s[i] ? -1.0 : 1.0;

      m0 = {0.0, 0.5};
      m1 = {1.0, 0.5};
      m2 = {0.5, 0.0};
      m3 = {0.5, 1.0};

      n0 = {-1.0,  0.0};
      n1 = { 1.0,  0.0};
      n2 = { 0.0, -1.0};
      n3 = { 0.0,  1.0};
    }

    template<typename F, typename C>
    void interpolate (const F& ff, std::vector<C>& out) const
    {
      // f gives v*outer normal at a point on the edge!
      auto&& f = Impl::makeFunctionWithCallOperator<typename LB::Traits::DomainType>(ff);

      out.resize(4);

      // Evaluate the normal components at the edge midpoints
      auto y = f(m0); out[0] = (y[0]*n0[0]+y[1]*n0[1])*sign_[0];
      y = f(m1); out[1] = (y[0]*n1[0]+y[1]*n1[1])*sign_[1];
      y = f(m2); out[2] = (y[0]*n2[0]+y[1]*n2[1])*sign_[2];
      y = f(m3); out[3] = (y[0]*n3[0]+y[1]*n3[1])*sign_[3];
    }

  private:
    std::array<typename LB::Traits::RangeFieldType,4> sign_;

    // The four edge midpoints of the reference quadrilateral
    typename LB::Traits::DomainType m0,m1,m2,m3;

    // The four edge normals of the reference quadrilateral
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
