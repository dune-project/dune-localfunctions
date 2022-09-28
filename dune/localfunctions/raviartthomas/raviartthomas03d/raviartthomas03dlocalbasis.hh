// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS_RAVIARTTHOMAS03D_RAVIARTTHOMAS03DLOCALBASIS_HH
#define DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS_RAVIARTTHOMAS03D_RAVIARTTHOMAS03DLOCALBASIS_HH

#include <numeric>

#include <dune/common/fmatrix.hh>

#include <dune/localfunctions/common/localbasis.hh>

namespace Dune
{
  /**@ingroup LocalBasisImplementation
         \brief Lowest order Raviart-Thomas shape functions on the reference tetrahedron.

         \tparam D Type to represent the field in the domain.
         \tparam R Type to represent the field in the range.

         \nosubgrouping
   */
  template<class D, class R>
  class RT03DLocalBasis
  {
  public:
    typedef LocalBasisTraits<D,3,Dune::FieldVector<D,3>,R,3,Dune::FieldVector<R,3>,
        Dune::FieldMatrix<R,3,3> > Traits;

    //! \brief Make set number s, where 0 <= s < 16
    RT03DLocalBasis (std::bitset<4> s = 0)
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
      auto c = std::sqrt(2.0);
      out[0] = {sign_[0]*c* in[0],       sign_[0]*c* in[1],       sign_[0]*c*(in[2]-D(1))};
      out[1] = {sign_[1]*c* in[0],       sign_[1]*c*(in[1]-D(1)), sign_[1]*c* in[2]      };
      out[2] = {sign_[2]*c*(in[0]-D(1)), sign_[2]*c* in[1],       sign_[2]*c* in[2]      };
      out[3] = {sign_[3]*c* in[0],       sign_[3]*c* in[1],       sign_[3]*c* in[2]      };
    }

    //! \brief Evaluate Jacobian of all shape functions
    inline void
    evaluateJacobian (const typename Traits::DomainType& in,                  // position
                      std::vector<typename Traits::JacobianType>& out) const  // return value
    {
      out.resize(4);
      for (int i=0; i<4; i++)
      {
        auto c = std::sqrt(2.0);
        out[i][0] = {c*sign_[i],         0,         0};
        out[i][1] = {         0,c*sign_[i],         0};
        out[i][2] = {         0,         0,c*sign_[i]};
      }
    }

    //! \brief Evaluate partial derivatives of all shape functions
    void partial (const std::array<unsigned int, 3>& order,
                  const typename Traits::DomainType& in,                 // position
                  std::vector<typename Traits::RangeType>& out) const    // return value
    {
      auto totalOrder = std::accumulate(order.begin(), order.end(), 0);
      if (totalOrder == 0) {
        evaluateFunction(in, out);
      } else if (totalOrder == 1) {
        auto const direction = std::distance(order.begin(), std::find(order.begin(), order.end(), 1));
        out.resize(size());

        for (int i=0; i<size(); i++)
        {
          out[i][direction] = sign_[i]* std::sqrt(2.0) ;
          out[i][(direction+1)%3] = 0;
          out[i][(direction+2)%3] = 0;
        }
      } else {
        out.resize(size());
        for (std::size_t i = 0; i < size(); ++i)
          for (std::size_t j = 0; j < 3; ++j)
            out[i][j] = 0;
      }

    }

    //! \brief Polynomial order of the shape functions
    unsigned int order () const
    {
      return 1;
    }

  private:

    // Signs of the face normals
    std::array<R,4> sign_;
  };
}
#endif
