// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_BREZZIDOUGLASFORTINMARINI_CUBE_LOCALBASIS_HH
#define DUNE_LOCALFUNCTIONS_BREZZIDOUGLASFORTINMARINI_CUBE_LOCALBASIS_HH

#include <algorithm>
#include <array>
#include <bitset>
#include <numeric>
#include <vector>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/math.hh>
#include <dune/common/rangeutilities.hh>
#include <dune/common/typetraits.hh>

#include <dune/localfunctions/common/localbasis.hh>

namespace Dune
{
  /**
   * \ingroup LocalBasisImplementation
   * \brief Brezzi-Douglas-Fortin-Marini shape functions on a reference cube.
   *
   * \tparam D      Type to represent the field in the domain.
   * \tparam R      Type to represent the field in the range.
   * \tparam dim    dimension of the reference element, must be >= 2.
   * \tparam order  order of the element, must be >= 1.
   *
   * \nosubgroup
   */
  template<class D, class R, unsigned int dim, unsigned int order>
  class BDFMCubeLocalBasis
  {
    static_assert( AlwaysFalse<D>::value,
                   "`BDFMCubeLocalBasis` not implemented for chosen `dim` and `order`." );
  };


#ifndef DOXYGEN
  template<class D, class R, unsigned int dim>
  class BDFMCubeLocalBasis<D, R, dim, 0>
  {
    static_assert( AlwaysFalse<D>::value,
                   "`BDFMCubeLocalBasis` not defined for order 0." );
  };
#endif // #ifndef DOXYGEN


  /**
   * \brief First order Brezzi-Douglas-Fortin-Marini shape functions on the refrence quadrialteral.
   *
   * \nosubgrouping
   */
  template<class D, class R >
  class BDFMCubeLocalBasis<D, R, 2, 1>
  {
    using DomainType    = FieldVector<D, 2>;
    using RangeType     = FieldVector<R, 2>;
    using JacobianType  = FieldMatrix<R, 2, 2>;

  public:
    using Traits = LocalBasisTraits<D, 2, DomainType, R, 2, RangeType, JacobianType>;

    //! \brief Standard constructor
    BDFMCubeLocalBasis ()
    {
      std::fill(s_.begin(), s_.end(), 1);
    }

    /**
     * \brief Make set number s, where 0<= s < 16
     *
     * \param s  Edge orientation indicator
     */
    BDFMCubeLocalBasis (std::bitset<4> s)
    {
      for (auto i : range(4))
        s_[i] = s[i] ? -1 : 1;
    }

    //! \brief number of shape functions
    unsigned int size () const { return 4; }

    /**
     * \brief Evaluate all shape functions
     *
     * \param in   Position
     * \param out  return value
     */
    inline void evaluateFunction (const DomainType& in, std::vector<RangeType>& out) const
    {
      out.resize(4);

      out[0] = {s_[0]*(in[0]-1), 0 };
      out[1] = {s_[1]*(in[0])  , 0 };
      out[2] = {0, s_[2]*(in[1]-1)};
      out[3] = {0, s_[3]*(in[1]) };
    }

    /**
     * \brief Evaluate Jacobian of all shape functions
     *
     * \param in   Position
     * \param out  return value
     */
    inline void evaluateJacobian (const DomainType& in, std::vector<JacobianType>& out) const
    {
      out.resize(4);

      out[0] = {{s_[0], 0}, {0, 0}};
      out[1] = {{s_[1], 0}, {0, 0}};
      out[2] = {{0, 0}, {0, s_[2]}};
      out[3] = {{0, 0}, {0, s_[3]}};
    }

    /**
     * \brief Evaluate all partial derivatives of all shape functions
     *
     * \param order  order the partial derivative
     * \param in     Position
     * \param out    return value
     */
    void partial (const std::array<unsigned int, 2>& order,
                  const DomainType& in,
                  std::vector<RangeType>& out) const
    {
      if (std::accumulate(order.begin(), order.end(), 0) == 0)
        evaluateFunction(in, out);
      else
        DUNE_THROW(NotImplemented, "Desired derivative order is not implemented");
    }

    //! \brief Polynomial order of the shape functions
    unsigned int order () const { return 1; }

  private:
    std::array<R, 4> s_;
  };


  /**
   * \brief Second order Brezzi-Douglas-Fortin-Marini shape functions on the refrence quadrialteral.
   *
   * \nosubgrouping
   */
  template<class D, class R >
  class BDFMCubeLocalBasis<D, R, 2, 2>
  {
    using DomainType    = FieldVector<D, 2>;
    using RangeType     = FieldVector<R, 2>;
    using JacobianType  = FieldMatrix<R, 2, 2>;

  public:
    using Traits = LocalBasisTraits<D, 2, DomainType, R, 2, RangeType, JacobianType>;

    //! \brief Standard constructor
    BDFMCubeLocalBasis ()
    {
      std::fill(s_.begin(), s_.end(), 1);
    }

    /**
     * \brief Make set number s, where 0<= s < 16
     *
     * \param s  Edge orientation indicator
     */
    BDFMCubeLocalBasis (std::bitset<4> s)
    {
      for (auto i : range(4))
        s_[i] = s[i] ? -1 : 1;
    }

    //! \brief number of shape functions
    unsigned int size () const { return 10; }

    /**
     * \brief Evaluate all shape functions
     *
     * \param in   Position
     * \param out  return value
     */
    inline void evaluateFunction (const DomainType& in, std::vector<RangeType>& out) const
    {
      out.resize(10);

      const auto& x = in[0];
      const auto& y = in[1];

      auto pre = 1-x;
      out[0] = {pre*s_[0]*(3*x-1), 0};
      out[1] = {    pre*3*(2*y-1), 0};

      pre = x;
      out[2] = {pre*s_[1]*(3*x-2), 0};
      out[3] = {pre*    3*(2*y-1), 0};

      pre = 1-y;
      out[4] = {0, pre*s_[2]*(3*y-1)};
      out[5] = {0,     pre*3*(2*x-1)};

      pre = y;
      out[6] = {0, pre*s_[3]*(3*y-2)};
      out[7] = {0,     pre*3*(2*x-1)};

      out[8] = {6*x*(1-x), 0};

      out[9] = {0, 6*y*(1-y)};
    }

    /**
     * \brief Evaluate Jacobian of all shape functions
     *
     * \param in   Position
     * \param out  return value
     */
    inline void evaluateJacobian (const DomainType& in, std::vector<JacobianType>& out) const
    {
      out.resize(10);

      const auto& x = in[0];
      const auto& y = in[1];

      out[0] = {{s_[0]*(4-6*x), 0}, {0, 0}};
      out[1] = {{3-6*y,     6-6*x}, {0, 0}};

      out[2] = {{s_[1]*(6*x-2), 0}, {0, 0}};
      out[3] = {{6*y-3,       6*x}, {0, 0}};

      out[4] = {{0, 0}, {0, s_[2]*(4-6*y)}};
      out[5] = {{0, 0}, {6-6*y,     3-6*x}};

      out[6] = {{0, 0}, {0, s_[3]*(6*y-2)}};
      out[7] = {{0, 0}, {6*y,       6*x-3}};

      out[8] = {{6-12*x, 0}, {0, 0}};

      out[9] = {{0, 0}, {0, 6-12*y}};
    }

    /**
     * \brief Evaluate all partial derivatives of all shape functions
     *
     * \param order  order the partial derivative
     * \param in     Position
     * \param out    return value
     */
    void partial (const std::array<unsigned int, 2>& order,
                  const DomainType& in,
                  std::vector<RangeType>& out) const
    {
      if (std::accumulate(order.begin(), order.end(), 0) == 0)
        evaluateFunction(in, out);
      else
        DUNE_THROW(NotImplemented, "Desired derivative order is not implemented");
    }

    //! \brief Polynomial order of the shape functions
    unsigned int order () const { return 2; }

  private:
    std::array<R, 4> s_;
  };


  /**
   * \brief Third order Brezzi-Douglas-Fortin-Marini shape functions on the refrence quadrialteral.
   *
   * \nosubgrouping
   */
  template<class D, class R >
  class BDFMCubeLocalBasis<D, R, 2, 3>
  {
    using DomainType    = FieldVector<D, 2>;
    using RangeType     = FieldVector<R, 2>;
    using JacobianType  = FieldMatrix<R, 2, 2>;

  public:
    using Traits = LocalBasisTraits<D, 2, DomainType, R, 2, RangeType, JacobianType>;

    //! \brief Standard constructor
    BDFMCubeLocalBasis ()
    {
      std::fill(s_.begin(), s_.end(), 1);
    }

    /**
     * \brief Make set number s, where 0<= s < 16
     *
     * \param s  Edge orientation indicator
     */
    BDFMCubeLocalBasis (std::bitset<4> s)
    {
      for (auto i : range(4))
        s_[i] = s[i] ? -1 : 1;
    }

    //! \brief number of shape functions
    unsigned int size () const { return 18; }

    /**
     * \brief Evaluate all shape functions
     *
     * \param in   Position
     * \param out  return value
     */
    inline void evaluateFunction (const DomainType& in, std::vector<RangeType>& out) const
    {
      out.resize(18);

      const auto& x = in[0];
      const auto& y = in[1];

      auto pre = 1-x;
      out[0]  = {pre*s_[0]*-1*(10*x*x-8*x+1), 0};
      out[1]  = {pre*      3*(1-3*x)*(2*y-1), 0};
      out[2]  = {pre* s_[0]*-5*(6*y*y-6*y+1), 0};

      pre = x;
      out[3]  = {pre*s_[1]*(10*x*x-12*x+3), 0};
      out[4]  = {pre*    3*(3*x-2)*(2*y-1), 0};
      out[5]  = {pre*s_[1]*5*(6*y*y-6*y+1), 0};

      pre = 1-y;
      out[6]  = {0, pre*s_[2]*-1*(10*y*y-8*y+1)};
      out[7]  = {0, pre*      3*(1-3*y)*(2*x-1)};
      out[8]  = {0, pre*s_[2]*-5*( 6*x*x-6*x+1)};

      pre = y;
      out[9]  = {0, pre*s_[3]*(10*y*y-12*y+3)};
      out[10] = {0, pre*    3*(3*y-2)*(2*x-1)};
      out[11] = {0, pre*s_[3]*5*(6*x*x-6*x+1)};

      pre = x*(1-x);
      out[12] = {pre* 6        , 0};
      out[14] = {pre*30*(2*x-1), 0};
      out[16] = {pre*18*(2*y-1), 0};

      pre = y*(1-y);
      out[13] = {0, pre* 6        };
      out[15] = {0, pre*18*(2*x-1)};
      out[17] = {0, pre*30*(2*y-1)};
    }

    /**
     * \brief Evaluate Jacobian of all shape functions
     *
     * \param in   Position
     * \param out  return value
     */
    inline void evaluateJacobian (const DomainType& in, std::vector<JacobianType>& out) const
    {
      out.resize(18);

      const auto& x = in[0];
      const auto& y = in[1];

      out[0]  = {{s_[0]*(30*x*x-36*x+9),                      0}, {0, 0}};
      out[1]  = {{  6*(6*x*y-3*x-4*y+2),        6*(3*x*x-4*x+1)}, {0, 0}};
      out[2]  = {{s_[0]*5*(6*y*y-6*y+1), s_[0]*30*(x-1)*(2*y-1)}, {0, 0}};

      out[3]  = {{  s_[1]*30*x*x-24*x+3,                  0}, {0, 0}};
      out[4]  = {{    6*(3*x-1)*(2*y-1),        6*x*(3*x-2)}, {0, 0}};
      out[5]  = {{s_[1]*5*(6*y*y-6*y+1), s_[1]*30*x*(2*y-1)}, {0, 0}};

      out[6]  = {{0, 0}, {                     0, s_[2]*(30*y*y-36*y+9)}};
      out[7]  = {{0, 0}, {       6*(3*y*y-4*y+1),   6*(6*y*x-3*y-4*x+2)}};
      out[8]  = {{0, 0}, {s_[2]*30*(y-1)*(2*x-1), s_[2]*5*(6*x*x-6*x+1)}};

      out[9]  = {{0, 0}, {                 0,   s_[3]*30*y*y-24*y+3}};
      out[10] = {{0, 0}, {       6*y*(3*y-2),     6*(3*y-1)*(2*x-1)}};
      out[11] = {{0, 0}, {s_[3]*30*y*(2*x-1), s_[3]*5*(6*x*x-6*x+1)}};

      out[12] = {{         -6*(2*x-1),          0}, {0, 0}};
      out[14] = {{  -30*(6*x*x-6*x+1),          0}, {0, 0}};
      out[16] = {{-18*(2*x-1)*(2*y-1), 36*x*(1-x)}, {0, 0}};

      out[13] = {{0, 0}, {         0,          -6*(2*y-1)}};
      out[15] = {{0, 0}, {36*y*(1-y), -18*(2*y-1)*(2*x-1)}};
      out[17] = {{0, 0}, {         0,   -30*(6*y*y-6*y+1)}};
    }

    /**
     * \brief Evaluate all partial derivatives of all shape functions
     *
     * \param order  order the partial derivative
     * \param in     Position
     * \param out    return value
     */
    void partial (const std::array<unsigned int, 2>& order,
                  const DomainType& in,
                  std::vector<RangeType>& out) const
    {
      if (std::accumulate(order.begin(), order.end(), 0) == 0)
        evaluateFunction(in, out);
      else
        DUNE_THROW(NotImplemented, "Desired derivative order is not implemented");
    }

    //! \brief Polynomial order of the shape functions
    unsigned int order () const { return 3; }

  private:
    std::array<R, 4> s_;
  };

} // namespace Dune

#endif // #ifndef DUNE_LOCALFUNCTIONS_BREZZIDOUGLASFORTINMARINI_CUBE_LOCALBASIS_HH
