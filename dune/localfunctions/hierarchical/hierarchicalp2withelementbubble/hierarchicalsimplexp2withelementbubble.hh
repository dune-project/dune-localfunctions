// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_HIERARCHICAL_SIMPLEX_P2_WITH_ELEMENT_BUBBLE_LOCALBASIS_HH
#define DUNE_HIERARCHICAL_SIMPLEX_P2_WITH_ELEMENT_BUBBLE_LOCALBASIS_HH

/** \file
    \brief Hierarchical p2 shape functions for the simplex
 */

#include <numeric>
#include <vector>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <dune/localfunctions/common/localbasis.hh>
#include <dune/localfunctions/common/localkey.hh>
#include <dune/localfunctions/common/localinterpolation.hh>

namespace Dune
{
  template<class D, class R, int dim>
  class HierarchicalSimplexP2WithElementBubbleLocalBasis
  {
  public:
    HierarchicalSimplexP2WithElementBubbleLocalBasis()
    {
      DUNE_THROW(Dune::NotImplemented,"HierarchicalSimplexP2LocalBasis not implemented for dim > 3.");
    }
  };

  /**@ingroup LocalBasisImplementation
     \brief Hierarchical P2 basis in 1d.

     The shape functions are associated to the following points:

     f_0 ~ (0.0)   // linear function
     f_1 ~ (1.0)   // linear function
     f_2 ~ (0.5)   // quadratic bubble

     \tparam D Type to represent the field in the domain.
     \tparam R Type to represent the field in the range.

     \nosubgrouping
   */
  template<class D, class R>
  class HierarchicalSimplexP2WithElementBubbleLocalBasis<D,R,1>
  {
  public:
    //! \brief export type traits for function signature
    typedef LocalBasisTraits<D,1,Dune::FieldVector<D,1>,R,1,Dune::FieldVector<R,1>,
        Dune::FieldMatrix<R,1,1> > Traits;

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

      out[0] = 1-in[0];
      out[1] = in[0];
      out[2] = 1-4*(in[0]-0.5)*(in[0]-0.5);
    }

    //! \brief Evaluate Jacobian of all shape functions
    inline void
    evaluateJacobian (const typename Traits::DomainType& in,         // position
                      std::vector<typename Traits::JacobianType>& out) const      // return value
    {
      out.resize(3);

      out[0][0][0] = -1;
      out[1][0][0] =  1;
      out[2][0][0] = 4-8*in[0];
    }

    //! \brief Evaluate partial derivatives of all shape functions
    void partial (const std::array<unsigned int, 1>& order,
                  const typename Traits::DomainType& in,         // position
                  std::vector<typename Traits::RangeType>& out) const      // return value
    {
      auto totalOrder = order[0];
      if (totalOrder == 0) {
        evaluateFunction(in, out);
      } else if (totalOrder == 1) {
        out.resize(size());
        out[0] = -1;
        out[1] =  1;
        out[2] = 4-8*in[0];
      } else if (totalOrder == 2) {
        out.resize(size());
        out[0] = 0;
        out[1] = 0;
        out[2] =-8;
      } else {
        out.resize(size());
        out[0] = out[1] = out[2] = 0;
      }
    }

    /** \brief Polynomial order of the shape functions  (2, in this case)
     */
    unsigned int order () const
    {
      return 2;
    }

  };

  /**@ingroup LocalBasisImplementation
     \brief Hierarchical P2 basis in 1d.

     The shape functions are associated to the following points:

     The functions are associated to points by:

     f_0 ~ (0.0, 0.0)
     f_1 ~ (0.5, 0.0)
     f_2 ~ (1.0, 0.0)
     f_3 ~ (0.0, 0.5)
     f_4 ~ (0.5, 0.5)
     f_5 ~ (0.0, 1.0)
     f_6 ~ (1/3, 1/3)

     \tparam D Type to represent the field in the domain.
     \tparam R Type to represent the field in the range.

     \nosubgrouping
   */
  template<class D, class R>
  class HierarchicalSimplexP2WithElementBubbleLocalBasis<D,R,2>
  {
  public:
    //! \brief export type traits for function signature
    typedef LocalBasisTraits<D,2,Dune::FieldVector<D,2>,R,1,Dune::FieldVector<R,1>,
        Dune::FieldMatrix<R,1,2> > Traits;

    //! \brief number of shape functions
    unsigned int size () const
    {
      return 7;
    }

    //! \brief Evaluate all shape functions
    inline void evaluateFunction (const typename Traits::DomainType& in,
                                  std::vector<typename Traits::RangeType>& out) const
    {
      out.resize(7);

      out[0] = 1 - in[0] - in[1];
      out[1] = 4*in[0]*(1-in[0]-in[1]);
      out[2] = in[0];
      out[3] = 4*in[1]*(1-in[0]-in[1]);
      out[4] = 4*in[0]*in[1];
      out[5] = in[1];
      out[6] = 27*in[0]*in[1]*(1-in[0]-in[1]);

    }

    //! \brief Evaluate Jacobian of all shape functions
    inline void
    evaluateJacobian (const typename Traits::DomainType& in,         // position
                      std::vector<typename Traits::JacobianType>& out) const      // return value
    {
      out.resize(7);

      out[0][0][0] = -1;                    out[0][0][1] = -1;
      out[1][0][0] =  4-8*in[0]-4*in[1];    out[1][0][1] = -4*in[0];
      out[2][0][0] =  1;                    out[2][0][1] =  0;
      out[3][0][0] = -4*in[1];              out[3][0][1] =  4-4*in[0]-8*in[1];
      out[4][0][0] =  4*in[1];              out[4][0][1] =  4*in[0];
      out[5][0][0] =  0;                    out[5][0][1] =  1;

      // Cubic bubble
      out[6][0][0] = 27 * in[1] * (1 - 2*in[0] - in[1]);
      out[6][0][1] = 27 * in[0] * (1 - 2*in[1] - in[0]);

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
        out.resize(size());
        auto const direction = std::distance(order.begin(), std::find(order.begin(), order.end(), 1));

        switch (direction) {
        case 0:
          out[0] = -1;
          out[1] =  4-8*in[0]-4*in[1];
          out[2] =  1;
          out[3] = -4*in[1];
          out[4] =  4*in[1];
          out[5] =  0;
          out[6] = 27 * in[1] * (1 - 2*in[0] - in[1]);
          break;
        case 1:
          out[0] = -1;
          out[1] = -4*in[0];
          out[2] =  0;
          out[3] =  4-4*in[0]-8*in[1];
          out[4] =  4*in[0];
          out[5] =  1;
          out[6] = 27 * in[0] * (1 - 2*in[1] - in[0]);
          break;
        default:
          DUNE_THROW(RangeError, "Component out of range.");
        }
      } else {
        DUNE_THROW(NotImplemented, "Desired derivative order is not implemented");
      }
    }

    /** \brief Polynomial order of the shape functions  (3 in this case)
     */
    unsigned int order () const
    {
      return 3;
    }

  };

  /**@ingroup LocalBasisImplementation
     \brief Hierarchical P2 basis in 1d.

     The shape functions are associated to the following points:

     The functions are associated to points by:

     f_0 ~ (0.0, 0.0, 0.0)
     f_1 ~ (0.5, 0.0, 0.0)
     f_2 ~ (1.0, 0.0, 0.0)
     f_3 ~ (0.0, 0.5, 0.0)
     f_4 ~ (0.5, 0.5, 0.0)
     f_5 ~ (0.0, 1.0, 0.0)
     f_6 ~ (0.0, 0.0, 0.5)
     f_7 ~ (0.5, 0.0, 0.5)
     f_8 ~ (0.0, 0.5, 0.5)
     f_9 ~ (0.0, 0.0, 1.0)
     f_10 ~ (1/3, 1/3, 1/3)

     \tparam D Type to represent the field in the domain.
     \tparam R Type to represent the field in the range.

     \nosubgrouping
   */
  template<class D, class R>
  class HierarchicalSimplexP2WithElementBubbleLocalBasis<D,R,3>
  {
  public:
    //! \brief export type traits for function signature
    typedef LocalBasisTraits<D,3,Dune::FieldVector<D,3>,R,1,Dune::FieldVector<R,1>,
        Dune::FieldMatrix<R,1,3> > Traits;

    //! \brief number of shape functions
    unsigned int size () const
    {
      return 11;
    }

    //! \brief Evaluate all shape functions
    void evaluateFunction (const typename Traits::DomainType& in,
                           std::vector<typename Traits::RangeType>& out) const
    {
      out.resize(10);

      out[0] = 1 - in[0] - in[1] - in[2];
      out[1] = 4 * in[0] * (1 - in[0] - in[1] - in[2]);
      out[2] = in[0];
      out[3] = 4 * in[1] * (1 - in[0] - in[1] - in[2]);
      out[4] = 4 * in[0] * in[1];
      out[5] = in[1];
      out[6] = 4 * in[2] * (1 - in[0] - in[1] - in[2]);
      out[7] = 4 * in[0] * in[2];
      out[8] = 4 * in[1] * in[2];
      out[9] = in[2];

      // quartic element bubble
      out[10] = 81*in[0]*in[1]*in[2]*(1-in[0]-in[1]-in[2]);
    }

    //! \brief Evaluate Jacobian of all shape functions
    void evaluateJacobian (const typename Traits::DomainType& in,         // position
                           std::vector<typename Traits::JacobianType>& out) const      // return value
    {
      out.resize(10);

      out[0][0][0] = -1;                           out[0][0][1] = -1;                            out[0][0][2] = -1;
      out[1][0][0] =  4-8*in[0]-4*in[1]-4*in[2];   out[1][0][1] = -4*in[0];                      out[1][0][2] = -4*in[0];
      out[2][0][0] =  1;                           out[2][0][1] =  0;                            out[2][0][2] =  0;
      out[3][0][0] = -4*in[1];                     out[3][0][1] =  4-4*in[0]-8*in[1]-4*in[2];    out[3][0][2] = -4*in[1];
      out[4][0][0] =  4*in[1];                     out[4][0][1] =  4*in[0];                      out[4][0][2] =  0;
      out[5][0][0] =  0;                           out[5][0][1] =  1;                            out[5][0][2] =  0;
      out[6][0][0] = -4*in[2];                     out[6][0][1] = -4*in[2];                      out[6][0][2] =  4-4*in[0]-4*in[1]-8*in[2];
      out[7][0][0] =  4*in[2];                     out[7][0][1] =  0;                            out[7][0][2] =  4*in[0];
      out[8][0][0] =  0;                           out[8][0][1] =  4*in[2];                      out[8][0][2] =  4*in[1];
      out[9][0][0] =  0;                           out[9][0][1] =  0;                            out[9][0][2] =  1;

      out[10][0][0] = 81 * in[1] * in[2] * (1 - 2*in[0] -   in[1] -   in[2]);
      out[10][0][1] = 81 * in[0] * in[2] * (1 -   in[0] - 2*in[1] -   in[2]);
      out[10][0][2] = 81 * in[0] * in[1] * (1 -   in[0] -   in[1] - 2*in[2]);
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
          out[0] = -1;
          out[1] =  4-8*in[0]-4*in[1]-4*in[2];
          out[2] =  1;
          out[3] = -4*in[1];
          out[4] =  4*in[1];
          out[5] =  0;
          out[6] = -4*in[2];
          out[7] =  4*in[2];
          out[8] =  0;
          out[9] =  0;
          out[10] = 81 * in[1] * in[2] * (1 - 2*in[0] -   in[1] -   in[2]);
          break;
        case 1:
          out[0] = -1;
          out[1] = -4*in[0];
          out[2] =  0;
          out[3] =  4-4*in[0]-8*in[1]-4*in[2];
          out[4] =  4*in[0];
          out[5] =  1;
          out[6] = -4*in[2];
          out[7] =  0;
          out[8] =  4*in[2];
          out[9] =  0;
          out[10] = 81 * in[0] * in[2] * (1 -   in[0] - 2*in[1] -   in[2]);
          break;
        case 2:
          out[0] = -1;
          out[1] = -4*in[0];
          out[2] =  0;
          out[3] = -4*in[1];
          out[4] =  0;
          out[5] =  0;
          out[6] =  4-4*in[0]-4*in[1]-8*in[2];
          out[7] =  4*in[0];
          out[8] =  4*in[1];
          out[9] =  1;
          out[10] = 81 * in[0] * in[1] * (1 -   in[0] -   in[1] - 2*in[2]);
          break;
        default:
          DUNE_THROW(RangeError, "Component out of range.");
        }
      } else {
        DUNE_THROW(NotImplemented, "Desired derivative order is not implemented");
      }
    }

    /** \brief Polynomial order of the shape functions (4 in this case)
     */
    unsigned int order () const
    {
      return 4;
    }

  };


  /**@ingroup LocalBasisImplementation
     \brief The local finite element needed for the Zou-Kornhuber estimator for Signorini problems

     This shape function set consists of three parts:
     - Linear shape functions associated to the element vertices
     - Piecewise linear edge bubbles
     - A cubic element bubble

     Currently this element exists only for triangles!

     The functions are associated to points by:

     f_0 ~ (0.0, 0.0)
     f_1 ~ (1.0, 0.0)
     f_2 ~ (0.0, 1.0)
     f_3 ~ (0.5, 0.0)
     f_4 ~ (0.0, 0.5)
     f_5 ~ (0.5, 0.5)
     f_6 ~ (1/3, 1/3)

     \tparam D Type to represent the field in the domain.
     \tparam R Type to represent the field in the range.

     \nosubgrouping
   */
  template <int dim>
  class HierarchicalSimplexP2WithElementBubbleLocalCoefficients
  {
    // The binomial coefficient: dim+1 over 1
    static const int numVertices = dim+1;

    // The binomial coefficient: dim+1 over 2
    static const int numEdges    = (dim+1)*dim / 2;

  public:
    //! \brief Standard constructor
    HierarchicalSimplexP2WithElementBubbleLocalCoefficients ()
      : li(numVertices+numEdges + 1)
    {
      if (dim!=2)
        DUNE_THROW(NotImplemented, "only for 2d");

      li[0] = Dune::LocalKey(0,2,0);    // Vertex (0,0)
      li[1] = Dune::LocalKey(0,1,0);    // Edge   (0.5, 0)
      li[2] = Dune::LocalKey(1,2,0);    // Vertex (1,0)
      li[3] = Dune::LocalKey(1,1,0);    // Edge   (0, 0.5)
      li[4] = Dune::LocalKey(2,1,0);    // Edge   (0.5, 0.5)
      li[5] = Dune::LocalKey(2,2,0);    // Vertex (0,1)
      li[6] = Dune::LocalKey(0,0,0);    // Element (1/3, 1/3)
    }

    //! number of coefficients
    size_t size () const
    {
      return numVertices+numEdges + 1;
    }

    //! get i'th index
    const Dune::LocalKey& localKey (size_t i) const
    {
      return li[i];
    }

  private:
    std::vector<Dune::LocalKey> li;
  };

  template<class LB>
  class HierarchicalSimplexP2WithElementBubbleLocalInterpolation
  {
  public:

    //! \brief Local interpolation of a function
    template<typename F, typename C>
    void interpolate (const F& ff, std::vector<C>& out) const
    {
      typename LB::Traits::DomainType x;
      typename LB::Traits::RangeType y;

      out.resize(7);

      auto&& f = Impl::makeFunctionWithCallOperator<decltype(x)>(ff);

      // vertices
      x[0] = 0.0; x[1] = 0.0; out[0] = f(x);
      x[0] = 1.0; x[1] = 0.0; out[2] = f(x);
      x[0] = 0.0; x[1] = 1.0; out[5] = f(x);

      // edge bubbles
      x[0] = 0.5; x[1] = 0.0; y = f(x);
      out[1] = y - out[0]*(1-x[0]) - out[2]*x[0];

      x[0] = 0.0; x[1] = 0.5; y = f(x);
      out[3] = y - out[0]*(1-x[1]) - out[5]*x[1];

      x[0] = 0.5; x[1] = 0.5; y = f(x);
      out[4] = y - out[2]*x[0] - out[5]*x[1];

      // element bubble
      x[0] = 1.0/3; x[1] = 1.0/3; y = f(x);

      /** \todo Hack: extract the proper types */
      HierarchicalSimplexP2WithElementBubbleLocalBasis<double,double,2> shapeFunctions;
      std::vector<typename LB::Traits::RangeType> sfValues;
      shapeFunctions.evaluateFunction(x, sfValues);

      out[6] = y;
      for (int i=0; i<6; i++)
        out[6] -= out[i]*sfValues[i];

    }

  };


}
#endif
