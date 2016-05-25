// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS12DLOCALBASIS_HH
#define DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS12DLOCALBASIS_HH

#include <vector>

#include <dune/common/fmatrix.hh>

#include "../../common/localbasis.hh"

namespace Dune
{
  /**
   * @ingroup LocalBasisImplementation
   * \brief First order Raviart-Thomas shape functions on the reference triangle.
   *
   * \tparam D Type to represent the field in the domain.
   * \tparam R Type to represent the field in the range.
   *
   * \nosubgrouping
   */
  template<class D, class R>
  class RT12DLocalBasis
  {

  public:
    typedef LocalBasisTraits<D,2,Dune::FieldVector<D,2>,R,2,Dune::FieldVector<R,2>,
        Dune::FieldMatrix<R,2,2> > Traits;
    //! \brief Standard constructor
    RT12DLocalBasis ()
    {
      sign0 = sign1 = sign2 = 1.0;
    }

    /**
     * \brief Make set number s, where 0 <= s < 8
     *
     * \param s Edge orientation indicator
     */
    RT12DLocalBasis (unsigned int s)
    {
      sign0 = sign1 = sign2 = 1.0;
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
    }

    //! \brief number of shape functions
    constexpr unsigned int size () const
    {
      return 8;
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
      out.resize(8);
      out[0][0] = sign0*(in[0] - 4.0*in[0]*in[1]);
      out[0][1] = sign0*(-1.0 + 5.0*in[1] - 4.0*in[1]*in[1]);
      out[1][0] = sign1*(-1.0 + 5.0*in[0] - 4.0*in[0]*in[0]);
      out[1][1] = sign1*(in[1] - 4.0*in[0]*in[1]);
      out[2][0] = sign2*(-3.0*in[0] + 4.0*in[0]*in[0] + 4.0*in[1]*in[0]);
      out[2][1] = sign2*(-3.0*in[1] + 4.0*in[0]*in[1] + 4.0*in[1]*in[1]);
      out[3][0] = -5.0*in[0] + 8.0*in[0]*in[0] + 4.0*in[1]*in[0];
      out[3][1] = 3.0 - 6.0*in[0] - 7.0*in[1] + 8.0*in[0]*in[1] + 4.0*in[1]*in[1];
      out[4][0] = -3.0 + 7.0*in[0] + 6.0*in[1] - 4.0*in[0]*in[0] - 8.0*in[1]*in[0];
      out[4][1] = 5.0*in[1] - 4.0*in[0]*in[1] - 8.0*in[1]*in[1];
      out[5][0] = in[0] - 4.0*in[0]*in[0] + 4.0*in[1]*in[0];
      out[5][1] = -1.0*in[1] - 4.0*in[0]*in[1] + 4.0*in[1]*in[1];
      out[6][0] = 16.0*in[0] - 16.0*in[0]*in[0] - 8.0*in[1]*in[0];
      out[6][1] = 8.0*in[1] - 16.0*in[0]*in[1] - 8.0*in[1]*in[1];
      out[7][0] = 8.0*in[0] - 8.0*in[0]*in[0] - 16.0*in[1]*in[0];
      out[7][1] = 16.0*in[1] - 8.0*in[0]*in[1] - 16.0*in[1]*in[1];
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
      out.resize(8);

      out[0][0][0] = sign0*(1.0 - 4.0*in[1]);
      out[0][0][1] = sign0*(-4.0*in[0]);
      out[0][1][0] = 0.0;
      out[0][1][1] = sign0*(5.0 - 8.0*in[1]);

      out[1][0][0] = sign1*(5.0 - 8.0*in[0]);
      out[1][0][1] = 0.0;
      out[1][1][0] = sign1*(-4.0*in[1]);
      out[1][1][1] = sign1*(1.0 - 4.0*in[0]);

      out[2][0][0] = sign2*(-3.0 + 8.0*in[0] + 4.0*in[1]);
      out[2][0][1] = sign2*(4.0*in[0]);
      out[2][1][0] = sign2*(4.0*in[1]);
      out[2][1][1] = sign2*(-3.0 + 4.0*in[0] + 8.0*in[1]);

      out[3][0][0] = -5.0 + 16.0*in[0] + 4.0*in[1];
      out[3][0][1] = 4.0*in[0];
      out[3][1][0] = -6.0 + 8.0*in[1];
      out[3][1][1] = -7.0 + 8.0*in[0] + 8.0*in[1];

      out[4][0][0] = 7.0 - 8.0*in[0] - 8.0*in[1];
      out[4][0][1] = 6.0 - 8.0*in[0];
      out[4][1][0] = -4.0*in[1];
      out[4][1][1] = 5.0 - 4.0*in[0] - 16.0*in[1];

      out[5][0][0] = 1.0 - 8.0*in[0] + 4*in[1];
      out[5][0][1] = 4.0*in[0];
      out[5][1][0] = -4.0*in[1];
      out[5][1][1] = -1.0 - 4.0*in[0] + 8.0*in[1];

      out[6][0][0] = 16.0 - 32.0*in[0] - 8.0*in[1];
      out[6][0][1] = -8.0*in[0];
      out[6][1][0] = -16.0*in[1];
      out[6][1][1] = 8.0 - 16.0*in[0] - 16.0*in[1];

      out[7][0][0] = 8.0 - 16.0*in[0] - 16.0*in[1];
      out[7][0][1] = -16.0*in[0];
      out[7][1][0] = -8.0*in[1];
      out[7][1][1] = 16.0 - 8.0*in[0] - 32.0*in[1];
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
        // Calculate directions from order and call evaluate for the
        // specific totalOrder value, to calculate the derivatives.
        int dOrder = staticFindIf<1, Traits::diffOrder+1>([&](const auto i)
        {
          if (i == totalOrder) {
            std::array<int, i> directions;
            Impl::order2directions(order, directions);
            this->evaluate<i>(directions, in, out);
            return true; // terminate loop
          } else {
            return false;
          }
        });

        if (dOrder > Traits::diffOrder)
          DUNE_THROW(NotImplemented, "Desired derivative order is not implemented");
      }
    }

    template <std::size_t dOrder>
    inline void evaluate (const std::array<int, dOrder>& directions,
                          const typename Traits::DomainType& in,         // position
                          std::vector<typename Traits::RangeType>& out) const      // return value
    {
      if (dOrder == 0) {
        evaluateFunction(in, out);
      } else if (dOrder == 1) {
        out.resize(size());

        switch (directions[0]) {
          case 0:
            out[0][0] = sign0*(1.0 - 4.0*in[1]);
            out[0][1] = 0.0;
            out[1][0] = sign1*(5.0 - 8.0*in[0]);
            out[1][1] = sign1*(-4.0*in[1]);
            out[2][0] = sign2*(-3.0 + 8.0*in[0] + 4.0*in[1]);
            out[2][1] = sign2*(4.0*in[1]);
            out[3][0] = -5.0 + 16.0*in[0] + 4.0*in[1];
            out[3][1] = -6.0 + 8.0*in[1];
            out[4][0] = 7.0 - 8.0*in[0] - 8.0*in[1];
            out[4][1] = -4.0*in[1];
            out[5][0] = 1.0 - 8.0*in[0] + 4*in[1];
            out[5][1] = -4.0*in[1];
            out[6][0] = 16.0 - 32.0*in[0] - 8.0*in[1];
            out[6][1] = -16.0*in[1];
            out[7][0] = 8.0 - 16.0*in[0] - 16.0*in[1];
            out[7][1] = -8.0*in[1];
            break;
          case 1:
            out[2][1] = sign2*(-3.0 + 4.0*in[0] + 8.0*in[1]);
            out[2][0] = sign2*(4.0*in[0]);
            out[1][1] = sign1*(1.0 - 4.0*in[0]);
            out[1][0] = 0.0;
            out[0][0] = sign0*(-4.0*in[0]);
            out[0][1] = sign0*(5.0 - 8.0*in[1]);
            out[3][0] = 4.0*in[0];
            out[3][1] = -7.0 + 8.0*in[0] + 8.0*in[1];
            out[4][0] = 6.0 - 8.0*in[0];
            out[4][1] = 5.0 - 4.0*in[0] - 16.0*in[1];
            out[5][0] = 4.0*in[0];
            out[5][1] = -1.0 - 4.0*in[0] + 8.0*in[1];
            out[6][0] = -8.0*in[0];
            out[6][1] = 8.0 - 16.0*in[0] - 16.0*in[1];
            out[7][0] = -16.0*in[0];
            out[7][1] = 16.0 - 8.0*in[0] - 32.0*in[1];
            break;
          default:
            DUNE_THROW(RangeError, "Component out of range.");
        }
      } else {
        DUNE_THROW(NotImplemented, "Desired derivative order is not implemented");
      }
    }

    //! \brief Polynomial order of the shape functions
    unsigned int order () const
    {
      return 2;
    }

  private:
    R sign0, sign1, sign2;
  };
}
#endif // DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS12DLOCALBASIS_HH
