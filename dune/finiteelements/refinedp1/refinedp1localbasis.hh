// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_REFINED_P1_LOCALBASIS_HH
#define DUNE_REFINED_P1_LOCALBASIS_HH

/** \file
    \brief Linear Lagrange shape functions on a uniformly refined reference element
 */

#include "../common/localbasis.hh"

namespace Dune
{
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
  class RefinedP1LocalBasis :
    public C1LocalBasisInterface<
        C1LocalBasisTraits<D,2,Dune::FieldVector<D,2>,R,1,Dune::FieldVector<R,1>,
            Dune::FieldVector<Dune::FieldVector<R,2>,1> >
#ifndef DUNE_VIRTUAL_SHAPEFUNCTIONS
        ,RefinedP1LocalBasis<D,R>
#endif
        >
  {
  public:
    //! \brief export type traits for function signature
    typedef C1LocalBasisTraits<D,2,Dune::FieldVector<D,2>,R,1,Dune::FieldVector<R,1>,
        Dune::FieldVector<Dune::FieldVector<R,2>,1> > Traits;

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

      int subElement;
      typename Traits::DomainType local;
      getSubElement(in, subElement, local);

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
      getSubElement(in, subElement, local);

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

    /** \brief Polynomial order of the shape functions
        Doesn't really apply: these shape functions are only piecewise linear
     */
    unsigned int order () const
    {
      return 1;
    }

  private:
    /** \brief Get local coordinates in the subtriangle

       \param[in] global Coordinates in the reference triangle
       \param[out] subElement Which of the four subtriangles is global in?
       \param[out] local The local coordinates in the subtriangle
     */
    static void getSubElement(const typename Traits::DomainType& global,
                              int& subElement,
                              typename Traits::DomainType& local)
    {
      if (global[0] + global[1] <= 0.5) {
        subElement = 0;
        local[0] = 2*global[0];
        local[1] = 2*global[1];
        return;
      } else if (global[0] >= 0.5) {
        subElement = 1;
        local[0] = 2*global[0]-1;
        local[1] = 2*global[1];
        return;
      } else if (global[1] >= 0.5) {
        subElement = 2;
        local[0] = 2*global[0];
        local[1] = 2*global[1]-1;
        return;
      }

      subElement = 3;
      local[0] = -2 * global[0] + 1;
      local[1] = -2 * global[1] + 1;

    }

  };
}
#endif
