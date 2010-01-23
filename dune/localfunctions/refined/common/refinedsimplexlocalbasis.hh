// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_REFINED_SIMPLEX_LOCALBASIS_HH
#define DUNE_REFINED_SIMPLEX_LOCALBASIS_HH

/** \file
    \brief Contains a base class for a LocalBasis classes based on uniform refinement
 */

#include <dune/common/fvector.hh>
#include <dune/common/exceptions.hh>
#include <dune/localfunctions/common/localbasis.hh>

namespace Dune
{
  template<class D, int dim>
  class RefinedSimplexLocalBasis
  {
  protected:
    RefinedSimplexLocalBasis()
    {
      DUNE_THROW(Dune::NotImplemented,"RefinedSimplexLocalBasis not implemented for dim > 3.");
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

     \nosubgrouping
   */
  template<class D>
  class RefinedSimplexLocalBasis<D,1>
  {
  protected:

    /** \brief Protected default constructor so this class can only be instantiated as a base class. */
    RefinedSimplexLocalBasis() {}

    /** \brief Get local coordinates in the subtriangle

       \param[in] global Coordinates in the reference triangle
       \param[out] subElement Which of the four subtriangles is global in?
       \param[out] local The local coordinates in the subtriangle
     */
    static void getSubElement(const FieldVector<D,1>& global,
                              int& subElement,
                              FieldVector<D,1>& local)
    {
      if (global[0] <= 0.5) {
        subElement = 0;
        local[0] = 2.0 * global[0];
        return;
      }

      subElement = 1;
      local[0] = 2.0 * global[0] - 1.0;

    }

  };


  /**@ingroup LocalBasisImplementation
     \brief Uniformly refined constant shape functions on the triangle.

     This shape function set mimicks the P0 shape functions that you would get on
     a uniformly refined grid.  Hence these shape functions are only piecewise
     constant!

     Shape functions like these are necessary for hierarchical error estimators
     for certain nonlinear problems.

     The functions are associated with the simplices having the following centers:

     f_0 ~ (1/6, 1/6)
     f_1 ~ (4/6, 1/6)
     f_2 ~ (1/6, 4/6)
     f_3 ~ (2/6, 2/6)

     \tparam D Type to represent the field in the domain.

     \nosubgrouping
   */
  template<class D>
  class RefinedSimplexLocalBasis<D,2>
  {
  protected:

    /** \brief Protected default constructor so this class can only be instantiated as a base class. */
    RefinedSimplexLocalBasis() {}

    /** \brief Get local coordinates in the subtriangle.
     *
     * The triangles are ordered according to
     *
     * |\
     * |2\
     * |--\
     * |\3|\
     * |0\|1\
     * ------
     *
     * \param[in] global Coordinates in the reference triangle
     * \returns Number of the subtriangles containing in
     */
    static int getSubElement(const FieldVector<D,2>& global)
    {
      if (global[0] + global[1] <= 0.5)
        return 0;
      else if (global[0] >= 0.5)
        return 1;
      else if (global[1] >= 0.5)
        return 2;
      return 3;
    }

    /** \brief Get local coordinates in the subtriangle

       \param[in] global Coordinates in the reference triangle
       \param[out] subElement Which of the four subtriangles is global in?
       \param[out] local The local coordinates in the subtriangle
     */
    static void getSubElement(const FieldVector<D,2>& global,
                              int& subElement,
                              FieldVector<D,2>& local)
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

     \nosubgrouping
   */
  template<class D>
  class RefinedSimplexLocalBasis<D,3>
  {
  protected:

    /** \brief Protected default constructor so this class can only be instantiated as a base class. */
    RefinedSimplexLocalBasis() {}

    /** \brief Get local coordinates in the subsimplex

       \param[in] global Coordinates in the reference simplex
       \param[out] subElement Which of the subsimplice is global in?
       \param[out] local The local coordinates in the subsimplex
     */
    static void getSubElement(const FieldVector<D,3>& global,
                              int& subElement,
                              FieldVector<D,3>& local)
    {
      if (global[0] + global[1] + global[2] <= 0.5) {
        subElement = 0;
        local = global;
        local *= 2.0;
        return;
      } else if (global[0] >= 0.5) {
        subElement = 1;
        local = global;
        local[0] -= 0.5;
        local *= 2.0;
        return;
      } else if (global[1] >= 0.5) {
        subElement = 2;
        local = global;
        local[1] -= 0.5;
        local *= 2.0;
        return;
      } else if (global[2] >= 0.5) {
        subElement = 3;
        local = global;
        local[2] -= 0.5;
        local *= 2.0;
        return;
      } else if ((global[0] + global[1] <= 0.5)and (global[1] + global[2] <= 0.5)) {
        subElement = 4;
        local[0] = 2.0 * global[1];
        local[1] = 2.0 * (0.5 - global[0] - global[1]);
        local[2] = 2.0 * (-0.5 + global[0] + global[1] + global[2]);
        //              Dune::FieldMatrix<double,3,3> A(0.0);
        //              A[0][1] =  2.0;
        //              A[1][0] = -2.0;
        //              A[1][1] = -2.0;
        //              A[2][0] =  2.0;
        //              A[2][1] =  2.0;
        //              A[2][2] =  2.0;
        //              A.mv(global,local);
        //              local[1] += 1.0;
        //              local[2] -= 1.0;
        return;
      } else if ((global[0] + global[1] >= 0.5)and (global[1] + global[2] <= 0.5)) {
        subElement = 5;
        local[0] = 2.0 * (0.5 - global[0]);
        local[1] = 2.0 * (0.5 - global[1] - global[2]);
        local[2] = 2.0 * global[2];
        //              Dune::FieldMatrix<double,3,3> A(0.0);
        //              A[0][0] = -2.0;
        //              A[1][1] = -2.0;
        //              A[1][2] = -2.0;
        //              A[2][2] =  2.0;
        //              A.mv(global,local);
        //              local[0] += 1.0;
        //              local[1] += 1.0;
        return;
      } else if ((global[0] + global[1] <= 0.5)and (global[1] + global[2] >= 0.5)) {
        subElement = 6;
        local[0] = 2.0 * (0.5 - global[0] - global[1]);
        local[1] = 2.0 * global[0];
        local[2] = 2.0 * (-0.5 + global[1] + global[2]);
        //              Dune::FieldMatrix<double,3,3> A(0.0);
        //              A[0][0] = -2.0;
        //              A[0][1] = -2.0;
        //              A[1][0] =  2.0;
        //              A[2][1] =  2.0;
        //              A[2][2] =  2.0;
        //              A.mv(global,local);
        //              local[0] += 1.0;
        //              local[2] -= 1.0;
        return;
      } else if ((global[0] + global[1] >= 0.5)and (global[1] + global[2] >= 0.5)) {
        subElement = 7;
        local[0] = 2.0 * (-0.5 + global[1] + global[2]);
        local[1] = 2.0 * (0.5 - global[1]);
        local[2] = 2.0 * (-0.5 + global[0] + global[1]);
        //              Dune::FieldMatrix<double,3,3> A(0.0);
        //              A[0][1] =  2.0;
        //              A[0][2] =  2.0;
        //              A[1][1] = -2.0;
        //              A[2][0] =  2.0;
        //              A[2][1] =  2.0;
        //              A.mv(global,local);
        //              local[0] -= 1.0;
        //              local[1] += 1.0;
        //              local[2] -= 1.0;
        return;
      }

      DUNE_THROW(InvalidStateException, "no subelement defined");

    }

  };


}

#endif
