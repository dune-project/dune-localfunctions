// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_REFINED_P1_LOCALBASIS_HH
#define DUNE_REFINED_P1_LOCALBASIS_HH

/** \file
    \brief Linear Lagrange shape functions on a uniformly refined reference element
 */

#include <dune/common/fmatrix.hh>

#include <dune/localfunctions/refined/common/refinedsimplexlocalbasis.hh>

namespace Dune
{
  template<class D, class R, int dim>
  class RefinedP1LocalBasis
    : public RefinedSimplexLocalBasis<D,dim>
  {
  public:
    RefinedP1LocalBasis()
    {
      DUNE_THROW(Dune::NotImplemented,"RefinedP1LocalBasis not implemented for dim > 3.");
    }
  };

  /**@ingroup LocalBasisImplementation
     \brief Uniformly refined linear Lagrange shape functions in 1D.

     1D IMPLEMENTATION IS NOT TESTED (the LocalElement for 1D does not exist due to lack of P21D elements)

     This shape function set mimicks the P1 shape functions that you would get on
     a uniformly refined grid.  Hence these shape functions are only piecewise
     linear!  The data layout is identical to P2 shape functions.

     Shape functions like these are necessary for hierarchical error estimators
     for certain nonlinear problems.

     The functions are associated to points by:

     f_0 ~ (0.0)
     f_1 ~ (1.0)
     f_2 ~ (0.5)

     \tparam D Type to represent the field in the domain.
     \tparam R Type to represent the field in the range.

     \nosubgrouping
   */
  template<class D, class R>
  class RefinedP1LocalBasis<D,R,1>
    : public RefinedSimplexLocalBasis<D,1>
  {
  public:
    //! \brief export type traits for function signature
    typedef LocalBasisTraits<D,1,Dune::FieldVector<D,1>,R,1,Dune::FieldVector<R,1>,
        Dune::FieldMatrix<R,1,1>, DUNE_MAX_DIFF_ORDER > Traits;

    //! \brief number of shape functions
    constexpr std::size_t size () const
    {
      return 3;
    }

    //! \brief Evaluate all shape functions
    inline void evaluateFunction (const typename Traits::DomainType& in,
                                  std::vector<typename Traits::RangeType>& out) const
    {
      out.resize(3);

      int subElement;
      typename Traits::DomainType local;
      this->getSubElement(in, subElement, local);

      switch (subElement) {
      case 0 :

        out[0] = 1 - local[0];
        out[1] = local[0];
        out[2] = 0;
        break;

      case 1 :

        out[0] = 0;
        out[1] = 1 - local[0];
        out[2] = local[0];
        break;

      }

    }

    //! \brief Evaluate Jacobian of all shape functions
    inline void
    evaluateJacobian (const typename Traits::DomainType& in,         // position
                      std::vector<typename Traits::JacobianType>& out) const      // return value
    {
      out.resize(3);

      int subElement;
      typename Traits::DomainType local;
      this->getSubElement(in, subElement, local);

      switch (subElement) {
      case 0 :

        out[0][0][0] = -2;
        out[1][0][0] =  2;
        out[2][0][0] =  0;
        break;

      case 1 :

        out[0][0][0] =  0;
        out[1][0][0] = -2;
        out[2][0][0] =  2;
        break;

      }
    }

    //! \brief Evaluate partial derivatives of all shape functions
    inline void partial (const std::array<unsigned int, 1>& order,
                         const typename Traits::DomainType& in,         // position
                         std::vector<typename Traits::RangeType>& out) const      // return value
    {
      auto totalOrder = order[0];
      if (totalOrder == 0) {
        evaluateFunction(in, out);
      } else if (totalOrder == 1) {
        out.resize(3);

        int subElement;
        typename Traits::DomainType local;
        this->getSubElement(in, subElement, local);

        switch (subElement) {
          case 0:
            out[0] = -2;
            out[1] =  2;
            out[2] =  0;
            break;
          case 1:
            out[0] =  0;
            out[1] = -2;
            out[2] =  2;
            break;
        }
      } else {
        out.resize(3);
        out[0] = out[1] = out[2] = 0;
      }
    }

    //! \brief Evaluate partial derivatives of all shape functions
    template <std::size_t dOrder>
    inline void evaluate(const std::array<int, dOrder>& directions,
                         const typename Traits::DomainType& in,         // position
                         std::vector<typename Traits::RangeType>& out) const      // return value
    {
      if (dOrder == 0) {
        evaluateFunction(in, out);
      } else if (dOrder == 1) {
        out.resize(3);

        int subElement;
        typename Traits::DomainType local;
        this->getSubElement(in, subElement, local);

        switch (subElement) {
          case 0:
            out[0] = -2;
            out[1] =  2;
            out[2] =  0;
            break;
          case 1:
            out[0] =  0;
            out[1] = -2;
            out[2] =  2;
            break;
        }
      } else {
        out.resize(3);
        out[0] = out[1] = out[2] = 0;
      }
    }

    /** \brief Polynomial order of the shape functions
        Doesn't really apply: these shape functions are only piecewise linear
     */
    unsigned int order () const
    {
      return 1;
    }

  };

  /**@ingroup LocalBasisImplementation
     \brief Uniformly refined linear Lagrange shape functions on the triangle.

     This shape function set mimicks the P1 shape functions that you would get on
     a uniformly refined grid.  Hence these shape functions are only piecewise
     linear!  The data layout is identical to P2 shape functions.

     Shape functions like these are necessary for hierarchical error estimators
     for certain nonlinear problems.

     The functions are associated to points by:

     f_0 ~ (0.0, 0.0)
     f_1 ~ (0.5, 0.0)
     f_2 ~ (1.0, 0.0)
     f_3 ~ (0.0, 0.5)
     f_4 ~ (0.5, 0.5)
     f_5 ~ (0.0, 1.0)

     \tparam D Type to represent the field in the domain.
     \tparam R Type to represent the field in the range.

     \nosubgrouping
   */
  template<class D, class R>
  class RefinedP1LocalBasis<D,R,2>
    : public RefinedSimplexLocalBasis<D,2>
  {
  public:
    //! \brief export type traits for function signature
    typedef LocalBasisTraits<D,2,Dune::FieldVector<D,2>,R,1,Dune::FieldVector<R,1>,
        Dune::FieldMatrix<R,1,2>, 0 /*DUNE_MAX_DIFF_ORDER*/ > Traits;

    //! \brief number of shape functions
    constexpr std::size_t size () const
    {
      return 6;
    }

    //! \brief Evaluate all shape functions
    inline void evaluateFunction (const typename Traits::DomainType& in,
                                  std::vector<typename Traits::RangeType>& out) const
    {
      out.resize(6);

      int subElement;
      typename Traits::DomainType local;
      this->getSubElement(in, subElement, local);

      switch (subElement) {
      case 0 :

        out[0] = 1 - local[0] - local[1];
        out[1] = local[0];
        out[2] = 0;
        out[3] = local[1];
        out[4] = 0;
        out[5] = 0;
        break;

      case 1 :

        out[0] = 0;
        out[1] = 1 - local[0] - local[1];
        out[2] = local[0];
        out[3] = 0;
        out[4] = local[1];
        out[5] = 0;
        break;

      case 2 :

        out[0] = 0;
        out[1] = 0;
        out[2] = 0;
        out[3] = 1 - local[0] - local[1];
        out[4] = local[0];
        out[5] = local[1];
        break;
      case 3 :

        out[0] = 0;
        out[1] = local[1];
        out[2] = 0;
        out[3] = local[0];
        out[4] = 1 - local[0] - local[1];
        out[5] = 0;
      }

    }

    //! \brief Evaluate Jacobian of all shape functions
    inline void
    evaluateJacobian (const typename Traits::DomainType& in,         // position
                      std::vector<typename Traits::JacobianType>& out) const      // return value
    {
      out.resize(6);

      int subElement;
      typename Traits::DomainType local;
      this->getSubElement(in, subElement, local);

      switch (subElement) {
      case 0 :

        out[0][0][0] = -2;    out[0][0][1] = -2;
        out[1][0][0] =  2;    out[1][0][1] =  0;
        out[2][0][0] =  0;    out[2][0][1] =  0;
        out[3][0][0] =  0;    out[3][0][1] =  2;
        out[4][0][0] =  0;    out[4][0][1] =  0;
        out[5][0][0] =  0;    out[5][0][1] =  0;
        break;

      case 1 :

        out[0][0][0] =  0;    out[0][0][1] =  0;
        out[1][0][0] = -2;    out[1][0][1] = -2;
        out[2][0][0] =  2;    out[2][0][1] =  0;
        out[3][0][0] =  0;    out[3][0][1] =  0;
        out[4][0][0] =  0;    out[4][0][1] =  2;
        out[5][0][0] =  0;    out[5][0][1] =  0;
        break;

      case 2 :

        out[0][0][0] =  0;    out[0][0][1] =  0;
        out[1][0][0] =  0;    out[1][0][1] =  0;
        out[2][0][0] =  0;    out[2][0][1] =  0;
        out[3][0][0] = -2;    out[3][0][1] = -2;
        out[4][0][0] =  2;    out[4][0][1] =  0;
        out[5][0][0] =  0;    out[5][0][1] =  2;
        break;
      case 3 :

        out[0][0][0] =  0;    out[0][0][1] =  0;
        out[1][0][0] =  0;    out[1][0][1] = -2;
        out[2][0][0] =  0;    out[2][0][1] =  0;
        out[3][0][0] = -2;    out[3][0][1] =  0;
        out[4][0][0] =  2;    out[4][0][1] =  2;
        out[5][0][0] =  0;    out[5][0][1] =  0;
      }
    }

    //! \brief Evaluate partial derivatives of all shape functions
    inline void partial (const std::array<unsigned int, 2>& order,
                         const typename Traits::DomainType& in,         // position
                         std::vector<typename Traits::RangeType>& out) const      // return value
    {
      auto totalOrder = std::accumulate(order.begin(), order.end(), 0);
      if (totalOrder == 0) {
        evaluateFunction(in, out);
      } else if (totalOrder == 1) {
        DUNE_THROW(NotImplemented, "Desired derivative order is not implemented");
      } else {
        out.resize(size());
        for (std::size_t i = 0; i < size(); ++i)
          out[i] = 0;
      }
    }

    //! \brief Evaluate partial derivatives of all shape functions
    template <std::size_t dOrder>
    inline void evaluate(const std::array<int, dOrder>& directions,
                         const typename Traits::DomainType& in,         // position
                         std::vector<typename Traits::RangeType>& out) const      // return value
    {
      if (dOrder == 0) {
        evaluateFunction(in, out);
      } else if (dOrder == 1) {
        DUNE_THROW(NotImplemented, "Desired derivative order is not implemented");
      } else {
        out.resize(size());
        for (std::size_t i = 0; i < size(); ++i)
          out[i] = 0;
      }
    }

    /** \brief Polynomial order of the shape functions
        Doesn't really apply: these shape functions are only piecewise linear
     */
    unsigned int order () const
    {
      return 1;
    }

  };

  /**@ingroup LocalBasisImplementation
     \brief Uniformly refined linear Lagrange shape functions on the 3D-simplex (tetrahedron).

     This shape function set mimicks the P1 shape functions that you would get on
     a uniformly refined grid.  Hence these shape functions are only piecewise
     linear!  The data layout is identical to P2 shape functions.

     Shape functions like these are necessary for hierarchical error estimators
     for certain nonlinear problems.

     The functions are associated to points by:

     f_0 ~ (0.0, 0.0, 0.0)
     f_1 ~ (1.0, 0.0, 0.0)
     f_2 ~ (0.0, 1.0, 0.0)
     f_3 ~ (0.0, 0.0, 1.0)
     f_4 ~ (0.5, 0.0, 0.0)
     f_5 ~ (0.5, 0.5, 0.0)
     f_6 ~ (0.0, 0.5, 0.0)
     f_7 ~ (0.0, 0.0, 0.5)
     f_8 ~ (0.5, 0.0, 0.5)
     f_9 ~ (0.0, 0.5, 0.5)

     \tparam D Type to represent the field in the domain.
     \tparam R Type to represent the field in the range.

     \nosubgrouping
   */
  template<class D, class R>
  class RefinedP1LocalBasis<D,R,3>
    : public RefinedSimplexLocalBasis<D,3>
  {
  public:
    //! \brief export type traits for function signature
    typedef LocalBasisTraits<D,3,Dune::FieldVector<D,3>,R,1,Dune::FieldVector<R,1>,
        Dune::FieldMatrix<R,1,3>, 0 /*DUNE_MAX_DIFF_ORDER*/ > Traits;

    //! \brief number of shape functions
    constexpr std::size_t size () const
    {
      return 10;
    }

    //! \brief Evaluate all shape functions
    inline void evaluateFunction (const typename Traits::DomainType& in,
                                  std::vector<typename Traits::RangeType>& out) const
    {
      out.resize(10);

      int subElement;
      typename Traits::DomainType local;
      this->getSubElement(in, subElement, local);

      switch (subElement) {
      case 0 :

        out[0] = 1 - local[0] - local[1] - local[2];
        out[1] = local[0];
        out[2] = 0;
        out[3] = local[1];
        out[4] = 0;
        out[5] = 0;
        out[6] = local[2];
        out[7] = 0;
        out[8] = 0;
        out[9] = 0;
        break;

      case 1 :

        out[0] = 0;
        out[1] = 1 - local[0] - local[1] -local[2];
        out[2] = local[0];
        out[3] = 0;
        out[4] = local[1];
        out[5] = 0;
        out[6] = 0;
        out[7] = local[2];
        out[8] = 0;
        out[9] = 0;
        break;

      case 2 :

        out[0] = 0;
        out[1] = 0;
        out[2] = 0;
        out[3] = 1 - local[0] - local[1] -local[2];
        out[4] = local[0];
        out[5] = local[1];
        out[6] = 0;
        out[7] = 0;
        out[8] = local[2];
        out[9] = 0;
        break;

      case 3 :

        out[0] = 0;
        out[1] = 0;
        out[2] = 0;
        out[3] = 0;
        out[4] = 0;
        out[5] = 0;
        out[6] = 1 - local[0] - local[1] -local[2];
        out[7] = local[0];
        out[8] = local[1];
        out[9] = local[2];
        break;

      case 4 :

        out[0] = 0;
        out[1] = 1 - local[0] - local[1] -local[2];
        out[2] = 0;
        out[3] = local[0];
        out[4] = 0;
        out[5] = 0;
        out[6] = local[1];
        out[7] = local[2];
        out[8] = 0;
        out[9] = 0;
        break;

      case 5 :

        out[0] = 0;
        out[1] = local[1];
        out[2] = 0;
        out[3] = local[0];
        out[4] = 1 - local[0] - local[1] -local[2];
        out[5] = 0;
        out[6] = 0;
        out[7] = local[2];
        out[8] = 0;
        out[9] = 0;
        break;

      case 6 :

        out[0] = 0;
        out[1] = 0;
        out[2] = 0;
        out[3] = 1 - local[0] - local[1] -local[2];
        out[4] = 0;
        out[5] = 0;
        out[6] = local[0];
        out[7] = local[1];
        out[8] = local[2];
        out[9] = 0;
        break;

      case 7 :

        out[0] = 0;
        out[1] = 0;
        out[2] = 0;
        out[3] = 1 - local[0] - local[1] -local[2];
        out[4] = local[2];
        out[5] = 0;
        out[6] = 0;
        out[7] = local[1];
        out[8] = local[0];
        out[9] = 0;
        break;
      }

    }

    //! \brief Evaluate Jacobian of all shape functions
    inline void
    evaluateJacobian (const typename Traits::DomainType& in,         // position
                      std::vector<typename Traits::JacobianType>& out) const      // return value
    {
      out.resize(10);

      int subElement;
      typename Traits::DomainType local;
      this->getSubElement(in, subElement, local);

      switch (subElement) {
      case 0 :

        out[0][0][0] = -2;    out[0][0][1] = -2;    out[0][0][2] = -2;
        out[1][0][0] =  2;    out[1][0][1] =  0;    out[1][0][2] =  0;
        out[2][0][0] =  0;    out[2][0][1] =  0;    out[2][0][2] =  0;
        out[3][0][0] =  0;    out[3][0][1] =  2;    out[3][0][2] =  0;
        out[4][0][0] =  0;    out[4][0][1] =  0;    out[4][0][2] =  0;
        out[5][0][0] =  0;    out[5][0][1] =  0;    out[5][0][2] =  0;
        out[6][0][0] =  0;    out[6][0][1] =  0;    out[6][0][2] =  2;
        out[7][0][0] =  0;    out[7][0][1] =  0;    out[7][0][2] =  0;
        out[8][0][0] =  0;    out[8][0][1] =  0;    out[8][0][2] =  0;
        out[9][0][0] =  0;    out[9][0][1] =  0;    out[9][0][2] =  0;
        break;

      case 1 :

        out[0][0][0] =  0;    out[0][0][1] =  0;    out[0][0][2] =  0;
        out[1][0][0] = -2;    out[1][0][1] = -2;    out[1][0][2] = -2;
        out[2][0][0] =  2;    out[2][0][1] =  0;    out[2][0][2] =  0;
        out[3][0][0] =  0;    out[3][0][1] =  0;    out[3][0][2] =  0;
        out[4][0][0] =  0;    out[4][0][1] =  2;    out[4][0][2] =  0;
        out[5][0][0] =  0;    out[5][0][1] =  0;    out[5][0][2] =  0;
        out[6][0][0] =  0;    out[6][0][1] =  0;    out[6][0][2] =  0;
        out[7][0][0] =  0;    out[7][0][1] =  0;    out[7][0][2] =  2;
        out[8][0][0] =  0;    out[8][0][1] =  0;    out[8][0][2] =  0;
        out[9][0][0] =  0;    out[9][0][1] =  0;    out[9][0][2] =  0;
        break;

      case 2 :

        out[0][0][0] =  0;    out[0][0][1] =  0;    out[0][0][2] =  0;
        out[1][0][0] =  0;    out[1][0][1] =  0;    out[1][0][2] =  0;
        out[2][0][0] =  0;    out[2][0][1] =  0;    out[2][0][2] =  0;
        out[3][0][0] = -2;    out[3][0][1] = -2;    out[3][0][2] = -2;
        out[4][0][0] =  2;    out[4][0][1] =  0;    out[4][0][2] =  0;
        out[5][0][0] =  0;    out[5][0][1] =  2;    out[5][0][2] =  0;
        out[6][0][0] =  0;    out[6][0][1] =  0;    out[6][0][2] =  0;
        out[7][0][0] =  0;    out[7][0][1] =  0;    out[7][0][2] =  0;
        out[8][0][0] =  0;    out[8][0][1] =  0;    out[8][0][2] =  2;
        out[9][0][0] =  0;    out[9][0][1] =  0;    out[9][0][2] =  0;
        break;

      case 3 :

        out[0][0][0] =  0;    out[0][0][1] =  0;    out[0][0][2] =  0;
        out[1][0][0] =  0;    out[1][0][1] =  0;    out[1][0][2] =  0;
        out[2][0][0] =  0;    out[2][0][1] =  0;    out[2][0][2] =  0;
        out[3][0][0] =  0;    out[3][0][1] =  0;    out[3][0][2] =  0;
        out[4][0][0] =  0;    out[4][0][1] =  0;    out[4][0][2] =  0;
        out[5][0][0] =  0;    out[5][0][1] =  0;    out[5][0][2] =  0;
        out[6][0][0] = -2;    out[6][0][1] = -2;    out[6][0][2] = -2;
        out[7][0][0] =  2;    out[7][0][1] =  0;    out[7][0][2] =  0;
        out[8][0][0] =  0;    out[8][0][1] =  2;    out[8][0][2] =  0;
        out[9][0][0] =  0;    out[9][0][1] =  0;    out[9][0][2] =  2;
        break;

      case 4 :

        out[0][0][0] =  0;    out[0][0][1] =  0;    out[0][0][2] =  0;
        out[1][0][0] =  0;    out[1][0][1] = -2;    out[1][0][2] = -2;
        out[2][0][0] =  0;    out[2][0][1] =  0;    out[2][0][2] =  0;
        out[3][0][0] =  0;    out[3][0][1] =  2;    out[3][0][2] =  0;
        out[4][0][0] =  0;    out[4][0][1] =  0;    out[4][0][2] =  0;
        out[5][0][0] =  0;    out[5][0][1] =  0;    out[5][0][2] =  0;
        out[6][0][0] = -2;    out[6][0][1] = -2;    out[6][0][2] =  0;
        out[7][0][0] =  2;    out[7][0][1] =  2;    out[7][0][2] =  2;
        out[8][0][0] =  0;    out[8][0][1] =  0;    out[8][0][2] =  0;
        out[9][0][0] =  0;    out[9][0][1] =  0;    out[9][0][2] =  0;
        break;

      case 5 :

        out[0][0][0] =  0;    out[0][0][1] =  0;    out[0][0][2] =  0;
        out[1][0][0] =  0;    out[1][0][1] = -2;    out[1][0][2] = -2;
        out[2][0][0] =  0;    out[2][0][1] =  0;    out[2][0][2] =  0;
        out[3][0][0] = -2;    out[3][0][1] =  0;    out[3][0][2] =  0;
        out[4][0][0] =  2;    out[4][0][1] =  2;    out[4][0][2] =  0;
        out[5][0][0] =  0;    out[5][0][1] =  0;    out[5][0][2] =  0;
        out[6][0][0] =  0;    out[6][0][1] =  0;    out[6][0][2] =  0;
        out[7][0][0] =  0;    out[7][0][1] =  0;    out[7][0][2] =  2;
        out[8][0][0] =  0;    out[8][0][1] =  0;    out[8][0][2] =  0;
        out[9][0][0] =  0;    out[9][0][1] =  0;    out[9][0][2] =  0;
        break;

      case 6 :

        out[0][0][0] =  0;    out[0][0][1] =  0;    out[0][0][2] =  0;
        out[1][0][0] =  0;    out[1][0][1] =  0;    out[1][0][2] =  0;
        out[2][0][0] =  0;    out[2][0][1] =  0;    out[2][0][2] =  0;
        out[3][0][0] =  0;    out[3][0][1] =  0;    out[3][0][2] = -2;
        out[4][0][0] =  0;    out[4][0][1] =  0;    out[4][0][2] =  0;
        out[5][0][0] =  0;    out[5][0][1] =  0;    out[5][0][2] =  0;
        out[6][0][0] = -2;    out[6][0][1] = -2;    out[6][0][2] =  0;
        out[7][0][0] =  2;    out[7][0][1] =  0;    out[7][0][2] =  0;
        out[8][0][0] =  0;    out[8][0][1] =  2;    out[8][0][2] =  2;
        out[9][0][0] =  0;    out[9][0][1] =  0;    out[9][0][2] =  0;
        break;

      case 7 :

        out[0][0][0] =  0;    out[0][0][1] =  0;    out[0][0][2] =  0;
        out[1][0][0] =  0;    out[1][0][1] =  0;    out[1][0][2] =  0;
        out[2][0][0] =  0;    out[2][0][1] =  0;    out[2][0][2] =  0;
        out[3][0][0] = -2;    out[3][0][1] = -2;    out[3][0][2] = -2;
        out[4][0][0] =  2;    out[4][0][1] =  2;    out[4][0][2] =  0;
        out[5][0][0] =  0;    out[5][0][1] =  0;    out[5][0][2] =  0;
        out[6][0][0] =  0;    out[6][0][1] =  0;    out[6][0][2] =  0;
        out[7][0][0] =  0;    out[7][0][1] = -2;    out[7][0][2] =  0;
        out[8][0][0] =  0;    out[8][0][1] =  2;    out[8][0][2] =  2;
        out[9][0][0] =  0;    out[9][0][1] =  0;    out[9][0][2] =  0;
        break;
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
      } else if (totalOrder == 1) {
        DUNE_THROW(NotImplemented, "Desired derivative order is not implemented");
      } else {
        out.resize(size());
        for (std::size_t i = 0; i < size(); ++i)
          out[i] = 0;
      }
    }

    //! \brief Evaluate partial derivatives of all shape functions
    template <std::size_t dOrder>
    inline void evaluate(const std::array<int, dOrder>& directions,
                         const typename Traits::DomainType& in,         // position
                         std::vector<typename Traits::RangeType>& out) const      // return value
    {
      if (dOrder == 0) {
        evaluateFunction(in, out);
      } else if (dOrder == 1) {
        DUNE_THROW(NotImplemented, "Desired derivative order is not implemented");
      } else {
        out.resize(size());
        for (std::size_t i = 0; i < size(); ++i)
          out[i] = 0;
      }
    }


    /** \brief Polynomial order of the shape functions
     *  Doesn't really apply: these shape functions are only piecewise linear
     */
    unsigned int order () const
    {
      return 1;
    }

  };
}
#endif
