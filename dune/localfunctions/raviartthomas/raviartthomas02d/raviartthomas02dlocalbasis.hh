// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_RT0TRIANGLELOCALBASIS_HH
#define DUNE_RT0TRIANGLELOCALBASIS_HH

#include <numeric>

#include <dune/common/fmatrix.hh>

#include <dune/localfunctions/common/localbasis.hh>

namespace Dune
{
  /**@ingroup LocalBasisImplementation
         \brief Lowest order Raviart-Thomas shape functions on the reference triangle.

         \tparam D Type to represent the field in the domain.
         \tparam R Type to represent the field in the range.

         \nosubgrouping
   */
  template<class D, class R>
  class RT02DLocalBasis
  {
  public:
    typedef LocalBasisTraits<D,2,Dune::FieldVector<D,2>,R,2,Dune::FieldVector<R,2>,
        Dune::FieldMatrix<R,2,2> > Traits;

    //! \brief Standard constructor
    RT02DLocalBasis ()
    {
      sign0 = sign1 = sign2 = 1.0;
    }

    //! \brief Make set numer s, where 0<=s<8
    RT02DLocalBasis (unsigned int s)
    {
      sign0 = sign1 = sign2 = 1.0;
      if (s&1) sign0 = -1.0;
      if (s&2) sign1 = -1.0;
      if (s&4) sign2 = -1.0;
    }

    //! \brief number of shape functions
    unsigned int size () const
    {
      return 3;
    }

    //! \brief Evaluate all shape functions
    inline void evaluateFunction (const typename Traits::DomainType& in,
                                  std::vector<typename Traits::RangeType>& out) const
    {
      out.resize(3);
      out[0][0] = sign0*in[0];       out[0][1]=sign0*(in[1]-D(1));
      out[1][0] = sign1*(in[0]-D(1)); out[1][1]=sign1*in[1];
      out[2][0] = sign2*in[0];       out[2][1]=sign2*in[1];
    }

    //! \brief Evaluate Jacobian of all shape functions
    inline void
    evaluateJacobian (const typename Traits::DomainType& in,             // position
                      std::vector<typename Traits::JacobianType>& out) const                          // return value
    {
      out.resize(3);
      out[0][0][0] = sign0;       out[0][0][1] = 0;
      out[0][1][0] = 0;           out[0][1][1] = sign0;
      out[1][0][0] = sign1;       out[1][0][1] = 0;
      out[1][1][0] = 0;           out[1][1][1] = sign1;
      out[2][0][0] = sign2;       out[2][0][1] = 0;
      out[2][1][0] = 0;           out[2][1][1] = sign2;
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

        out[0][direction] = sign0;
        out[0][1-direction] = 0;
        out[1][direction] = sign1;
        out[1][1-direction] = 0;
        out[2][direction] = sign2;
        out[2][1-direction] = 0;
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
    R sign0, sign1, sign2;
  };
}
#endif
