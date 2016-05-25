// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PYRAMID_P2_LOCALBASIS_HH
#define DUNE_PYRAMID_P2_LOCALBASIS_HH

#include <dune/common/fmatrix.hh>

#include <dune/localfunctions/common/localbasis.hh>

namespace Dune
{
  /**@ingroup LocalBasisImplementation
     \brief Quadratic Lagrange shape functions on the pyramid.

     Defines the quadratic shape functions on pyramid.
     Taken from
     Liping Liu, Kevin B. Davies, Michal Krizek, and Li Guan:
     On Higher Order Pyramidal Finite Elements
     Adv. Appl. Math. Mech., Vol. 3, No. 2, pp. 131-140, 2011

     \tparam D Type to represent the field in the domain.
     \tparam R Type to represent the field in the range.

     \nosubgrouping
   */
  template<class D, class R>
  class PyramidP2LocalBasis
  {
  public:
    //! \brief export type traits for function signature
    typedef LocalBasisTraits<D,3,Dune::FieldVector<D,3>,R,1,Dune::FieldVector<R,1>,
        Dune::FieldMatrix<R,1,3>, 0 > Traits;

    //! \brief number of shape functions
    constexpr unsigned int size () const
    {
      return 14;
    }

    //! \brief Evaluate all shape functions
    inline void evaluateFunction (const typename Traits::DomainType& in,
                                  std::vector<typename Traits::RangeType>& out) const
    {
      out.resize(14);

      // transform to reference element with base [-1,1]^2
      const R x = 2.0*in[0] + in[2] - 1.0;
      const R y = 2.0*in[1] + in[2] - 1.0;
      const R z = in[2];

      if (x > y)
      {
        // vertices
        out[0] = 0.25*(x + z)*(x + z - 1)*(y - z - 1)*(y - z);
        out[1] = -0.25*(x + z)*(y - z)*((x + z + 1)*(-y + z + 1) - 4*z) - z*(x - y);
        out[2] = 0.25*(x + z)*(y - z)*(y - z + 1)*(x + z - 1);
        out[3] = 0.25*(y - z)*(x + z)*(y - z + 1)*(x + z + 1);
        out[4] = z*(2*z - 1);

        // lower edges
        out[5] = -0.5*(y - z + 1)*(x + z - 1)*((y - 1)*(x + 1) + z*(x - y + z + 1));
        out[6] = -0.5*(y - z + 1)*(((x + z + 1)*(y - 1)*x - z) + z*(2*y + 1));
        out[7] = -0.5*(x + z - 1)*(((y - z - 1)*(x + 1)*y - z) + z*(2*x + 1));
        out[8] = -0.5*(y - z + 1)*(x + z - 1)*(x + 1)*y;

        // upper edges
        out[9] = z*(x + z - 1)*(y - z - 1);
        out[10] = -z*((x + z + 1)*(y - z - 1) + 4*z);
        out[11] = -z*(y - z + 1)*(x + z - 1);
        out[12] = z*(y - z + 1)*(x + z + 1);

        // base face
        out[13] = (y - z + 1)*(x + z - 1)*((y - 1)*(x + 1) + z*(x - y + z + 1));
      }
      else
      {
        // vertices
        out[0] = 0.25*(y + z)*(y + z - 1)*(x - z - 1)*(x - z);
        out[1] = -0.25*(x - z)*(y + z)*(x - z + 1)*(-y - z + 1);
        out[2] = 0.25*(x - z)*(y + z)*((x - z - 1)*(y + z + 1) + 4*z) + z*(x - y);
        out[3] = 0.25*(y + z)*(x - z)*(x - z + 1)*(y + z + 1);
        out[4] = z*(2*z - 1);

        // lower edges
        out[5] = -0.5*(y + z - 1)*(((x - z - 1)*(y + 1)*x - z) + z*(2*y + 1));
        out[6] = -0.5*(x - z + 1)*(y + z - 1)*(y + 1)*x;
        out[7] = -0.5*(x - z + 1)*(y + z - 1)*(x - 1)*y;
        out[8] = -0.5*(x - z + 1)*(((y + z + 1)*(x - 1)*y - z) + z*(2*x + 1));

        // upper edges
        out[9] = z*(y + z - 1)*(x - z - 1);
        out[10] = -z*(x - z + 1)*(y + z - 1);
        out[11] = -z*((y + z + 1)*(x - z - 1) + 4*z);
        out[12] = z*(x - z + 1)*(y + z + 1);

        // base face
        out[13] = (x - z + 1)*(y + z - 1)*((y + 1)*(x - 1) - z*(x - y - z - 1));
      }
    }

    //! \brief Evaluate Jacobian of all shape functions
    inline void evaluateJacobian (const typename Traits::DomainType& in,         // position
                                  std::vector<typename Traits::JacobianType>& out) const      // return value
    {
      out.resize(14);

      // transform to reference element with base [-1,1]^2
      const R x = 2.0*in[0] + in[2] - 1.0;
      const R y = 2.0*in[1] + in[2] - 1.0;
      const R z = in[2];

      // transformation of the gradient leads to a multiplication
      // with the Jacobian [2 0 0; 0 2 0; 1 1 1]
      if (x > y)
      {
        // vertices
        out[0][0][0] = 0.5*(y - z - 1)*(y - z)*(2*x + 2*z - 1);
        out[0][0][1] = 0.5*(x + z)*(x + z - 1)*(2*y - 2*z - 1);
        out[0][0][2] = 0.5*(out[0][0][0] + out[0][0][1])
                       + 0.25*((2*x + 2*z - 1)*(y - z - 1)*(y - z)
                               + (x + z)*(x + z - 1)*(-2*y + 2*z + 1));

        out[1][0][0] = 2*(-0.25*((y - z)*((x + z + 1)*(-y + z + 1) - 4*z)
                                 + (x + z)*(y - z)*(-y + z + 1)) - z);
        out[1][0][1] = 2*(-0.25*((x + z)*((x + z + 1)*(-y + z + 1) - 4*z)
                                 + (x + z)*(y - z)*(-(x + z + 1))) + z);
        out[1][0][2] = 0.5*(out[1][0][0] + out[1][0][1])
                       - 0.25*((y - z)*((x + z + 1)*(-y + z + 1) - 4*z)
                               - (x + z)*((x + z + 1)*(-y + z + 1) - 4*z)
                               + (x + z)*(y - z)*(x - y + 2*z - 2))
                       - (x - y);

        out[2][0][0] = 0.5*(y - z)*(y - z + 1)*(2*x + 2*z - 1);
        out[2][0][1] = 0.5*(x + z)*(2*y - 2*z + 1)*(x + z - 1);
        out[2][0][2] = 0.5*(out[2][0][0] + out[2][0][1])
                       + 0.25*((y - x - 2*z)*(y - z + 1)*(x + z - 1)
                               + (x + z)*(y - z)*(y - x - 2*z + 2));

        out[3][0][0] = 0.5*(y - z)*(2*x + 2*z + 1)*(y - z + 1);
        out[3][0][1] = 0.5*(2*y - 2*z + 1)*(x + z)*(x + z + 1);
        out[3][0][2] = 0.5*(out[3][0][0] + out[3][0][1])
                       + 0.25*((y - x - 2*z)*(y - z + 1)*(x + z + 1)
                               + (y - z)*(x + z)*(y - x - 2*z));

        out[4][0][0] = 0;
        out[4][0][1] = 0;
        out[4][0][2] = 4*z - 1;

        // lower edges
        out[5][0][0] = -((y - z + 1)*((y - 1)*(x + 1) + z*(x - y + z + 1))
                         + (y - z + 1)*(x + z - 1)*((y - 1) + z));
        out[5][0][1] = -((x + z - 1)*((y - 1)*(x + 1) + z*(x - y + z + 1))
                         + (y - z + 1)*(x + z - 1)*((x + 1) - z));
        out[5][0][2] = 0.5*(out[5][0][0] + out[5][0][1])
                       - 0.5*((-x + y - 2*z + 2)*((y - 1)*(x + 1) + z*(x - y + z + 1))
                              + (y - z + 1)*(x + z - 1)*(x - y + 2*z + 1));

        out[6][0][0] = -(y - z + 1)*(2*x + z + 1)*(y - 1);
        out[6][0][1] = -(((x + z + 1)*(y - 1)*x - z) + z*(2*y + 1)
                         + (y - z + 1)*((x + z + 1)*x + 2*z));
        out[6][0][2] = 0.5*(out[6][0][0] + out[6][0][1])
                       - 0.5*(-(((x + z + 1)*(y - 1)*x - z) + z*(2*y + 1))
                              + (y - z + 1)*(((y - 1)*x - 1) + 2*y + 1));

        out[7][0][0] = -(((y - z - 1)*(x + 1)*y - z) + z*(2*x + 1)
                         + (x + z - 1)*((y - z - 1)*y + 2*z));
        out[7][0][1] = -(x + z - 1)*(2*y - z - 1)*(x + 1);
        out[7][0][2] = 0.5*(out[7][0][0] + out[7][0][1])
                       - 0.5*(((y - z - 1)*(x + 1)*y - z) + z*(2*x + 1)
                              + (x + z - 1)*((-(x + 1)*y - 1) + 2*x + 1));

        out[8][0][0] = -(y - z + 1)*(2*x + z)*y;
        out[8][0][1] = -(2*y - z + 1)*(x + z - 1)*(x + 1);
        out[8][0][2] = 0.5*(out[8][0][0] + out[8][0][1])
                       - 0.5*(-x + y - 2*z + 2)*(x + 1)*y;

        // upper edges
        out[9][0][0] = 2*z*(y - z - 1);
        out[9][0][1] = 2*z*(x + z - 1);
        out[9][0][2] = 0.5*(out[9][0][0] + out[9][0][1])
                       + (x + z - 1)*(y - z - 1) + z*(-x + y - 2*z);

        out[10][0][0] = -2*z*(y - z - 1);
        out[10][0][1] = -2*z*(x + z + 1);
        out[10][0][2] = 0.5*(out[10][0][0] + out[10][0][1])
                        - ((x + z + 1)*(y - z - 1) + 4*z)
                        - z*(-x + y - 2*z + 2);

        out[11][0][0] = -2*z*(y - z + 1);
        out[11][0][1] = -2*z*(x + z - 1);
        out[11][0][2] = 0.5*(out[11][0][0] + out[11][0][1])
                        - (y - z + 1)*(x + z - 1) - z*(-x + y - 2*z + 2);

        out[12][0][0] = 2*z*(y - z + 1);
        out[12][0][1] = 2*z*(x + z + 1);
        out[12][0][2] = 0.5*(out[12][0][0] + out[12][0][1])
                        + (y - z + 1)*(x + z + 1) + z*(-x + y - 2*z);

        // base face
        out[13][0][0] = 2*((y - z + 1)*((y - 1)*(x + 1) + z*(x - y + z + 1))
                           + (y - z + 1)*(x + z - 1)*(y - 1 + z));
        out[13][0][1] = 2*((x + z - 1)*((y - 1)*(x + 1) + z*(x - y + z + 1))
                           + (y - z + 1)*(x + z - 1)*(x + 1 - z));
        out[13][0][2] = 0.5*(out[13][0][0] + out[13][0][1])
                        + ((-x + y - 2*z + 2)*((y - 1)*(x + 1) + z*(x - y + z + 1))
                           + (y - z + 1)*(x + z - 1)*(x - y + 2*z + 1));
      }
      else
      {
        // vertices
        out[0][0][0] = 0.5*(y + z)*(y + z - 1)*(2*x - 2*z - 1);
        out[0][0][1] = 0.5*(2*y + 2*z - 1)*(x - z - 1)*(x - z);
        out[0][0][2] = 0.5*(out[0][0][0] + out[0][0][1])
                       + 0.25*((2*y + 2*z - 1)*(x - z - 1)*(x - z)
                               + (y + z)*(y + z - 1)*(-2*x + 2*z + 1));

        out[1][0][0] = -0.5*(y + z)*(2*x - 2*z + 1)*(-y - z + 1);
        out[1][0][1] = -0.5*(x - z)*(x - z + 1)*(-2*y - 2*z + 1);
        out[1][0][2] = 0.5*(out[1][0][0] + out[1][0][1])
                       - 0.25*((x - y - 2*z)*(x - z + 1)*(-y - z + 1)
                               + (x - z)*(y + z)*(-x + y + 2*z - 2));

        out[2][0][0] = 0.5*((y + z)*((x - z - 1)*(y + z + 1) + 4*z)
                            + (x - z)*(y + z)*(y + z + 1) + 4*z);
        out[2][0][1] = 0.5*((x - z)*((x - z - 1)*(y + z + 1) + 4*z)
                            + (x - z)*(y + z)*(x - z - 1) - 4*z);
        out[2][0][2] = 0.5*(out[2][0][0] + out[2][0][1])
                       + 0.25*((x - y - 2*z)*((x - z - 1)*(y + z + 1) + 4*z)
                               + (x - z)*(y + z)*(x - y - 2*z + 2) + 4*(x - y));

        out[3][0][0] = 0.5*(y + z)*(2*x - 2*z + 1)*(y + z + 1);
        out[3][0][1] = 0.5*(x - z)*(x - z + 1)*(2*y + 2*z + 1);
        out[3][0][2] = 0.5*(out[3][0][0] + out[3][0][1])
                       + 0.25*((x - y - 2*z)*(x - z + 1)*(y + z + 1)
                               + (y + z)*(x - z)*(x - y - 2*z));

        out[4][0][0] = 0;
        out[4][0][1] = 0;
        out[4][0][2] = 4*z - 1;

        // lower edges
        out[5][0][0] = -(y + z - 1)*(2*x - z - 1)*(y + 1);
        out[5][0][1] = -(((x - z - 1)*(y + 1)*x - z) + z*(2*y + 1)
                         + (y + z - 1)*((x - z - 1)*x + 2*z));
        out[5][0][2] = 0.5*(out[5][0][0] + out[5][0][1])
                       - 0.5*((((x - z - 1)*(y + 1)*x - z) + z*(2*y + 1))
                              + (y + z - 1)*((-(y + 1)*x - 1) + 2*y + 1));

        out[6][0][0] = -(2*x - z + 1)*(y + z - 1)*(y + 1);
        out[6][0][1] = -(x - z + 1)*(2*y + z)*x;
        out[6][0][2] = 0.5*(out[6][0][0] + out[6][0][1])
                       - 0.5*(x - y - 2*z + 2)*(y + 1)*x;

        out[7][0][0] = -(2*x - z)*(y + z - 1)*y;
        out[7][0][1] = -(x - z + 1)*(2*y + z - 1)*(x - 1);
        out[7][0][2] = 0.5*(out[7][0][0] + out[7][0][1])
                       - 0.5*(x - y - 2*z + 2)*(x - 1)*y;

        out[8][0][0] = -(((y + z + 1)*(x - 1)*y - z) + z*(2*x + 1)
                         + (x - z + 1)*((y + z + 1)*y + 2*z));
        out[8][0][1] = -(x - z + 1)*(2*y + z + 1)*(x - 1);
        out[8][0][2] = 0.5*(out[8][0][0] + out[8][0][1])
                       - 0.5*(-(((y + z + 1)*(x - 1)*y - z) + z*(2*x + 1))
                              + (x - z + 1)*(((x - 1)*y - 1) + 2*x + 1));

        // upper edges
        out[9][0][0] = 2*z*(y + z - 1);
        out[9][0][1] = 2*z*(x - z - 1);
        out[9][0][2] = 0.5*(out[9][0][0] + out[9][0][1])
                       + (y + z - 1)*(x - z - 1) + z*(x - y - 2*z);

        out[10][0][0] = -2*z*(y + z - 1);
        out[10][0][1] = -2*z*(x - z + 1);
        out[10][0][2] = 0.5*(out[10][0][0] + out[10][0][1])
                        - (x - z + 1)*(y + z - 1) - z*(x - y - 2*z + 2);

        out[11][0][0] = -2*z*(y + z + 1);
        out[11][0][1] = -2*z*(x - z - 1);
        out[11][0][2] = 0.5*(out[11][0][0] + out[11][0][1])
                        - ((y + z + 1)*(x - z - 1) + 4*z) - z*(x - y - 2*z + 2);

        out[12][0][0] = 2*z*(y + z + 1);
        out[12][0][1] = 2*z*(x - z + 1);
        out[12][0][2] = 0.5*(out[12][0][0] + out[12][0][1])
                        + (x - z + 1)*(y + z + 1) + z*(x - y - 2*z);

        // base face
        out[13][0][0] = 2*((y + z - 1)*((y + 1)*(x - 1) - z*(x - y - z - 1))
                           + (x - z + 1)*(y + z - 1)*(y + 1 - z));
        out[13][0][1] = 2*((x - z + 1)*((y + 1)*(x - 1) - z*(x - y - z - 1))
                           + (x - z + 1)*(y + z - 1)*(x - 1 + z));
        out[13][0][2] = 0.5*(out[13][0][0] + out[13][0][1])
                        + (x - y - 2*z + 2)*((y + 1)*(x - 1) - z*(x - y - z - 1))
                        + (x - z + 1)*(y + z - 1)*(-(x - y - 2*z - 1));
      }
    }

    //! \brief Evaluate partial derivatives of all shape functions
    inline void partial (const std::array<unsigned int, 3>& order,
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
    constexpr unsigned int order () const
    {
      return 2;
    }
  };
}
#endif
