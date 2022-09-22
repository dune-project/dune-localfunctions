// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_REFINED_SIMPLEX_LOCALBASIS_HH
#define DUNE_REFINED_SIMPLEX_LOCALBASIS_HH

/** \file
    \brief Contains a base class for LocalBasis classes based on uniform refinement
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
     \brief Base class for LocalBasis classes based on uniform refinement in 1D; provides numbering and local coordinates of subelements

     \tparam D Type to represent the field in the domain.

     \nosubgrouping
   */
  template<class D>
  class RefinedSimplexLocalBasis<D,1>
  {
  protected:

    /** \brief Protected default constructor so this class can only be instantiated as a base class. */
    RefinedSimplexLocalBasis() {}

    /** \brief Get the number of the subelement containing a given point.
     *
     * The subelements are ordered according to
     *
     *     0       1
     * |-------:-------|
     *
     * \param[in] global Coordinates in the reference element
     * \returns Number of the subtriangle containing <tt>global</tt>
     */
    static int getSubElement(const FieldVector<D,1>& global)
    {
      if (global[0] <= 0.5)
        return 0;
      else if (global[0] <= 1.0)
        return 1;

      DUNE_THROW(InvalidStateException, "no subelement defined");
    }

    /** \brief Get local coordinates in the subelement

       \param[in] global Coordinates in the reference element
       \param[out] subElement Number of the subelement containing <tt>global</tt>
       \param[out] local The local coordinates in the subelement
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
     \brief Base class for LocalBasis classes based on uniform refinement in 2D; provides numbering and local coordinates of subelements

     Shape functions like these are necessary for hierarchical error estimators
     for certain nonlinear problems.

     \tparam D Type to represent the field in the domain.

     \nosubgrouping
   */
  template<class D>
  class RefinedSimplexLocalBasis<D,2>
  {
  protected:

    /** \brief Protected default constructor so this class can only be instantiated as a base class. */
    RefinedSimplexLocalBasis() {}

    /** \brief Get the number of the subtriangle containing a given point.
     *
     * The triangles are ordered according to
     * \verbatim
       |\
       |2\
       |--\
       |\3|\
       |0\|1\
       ------
       \endverbatim
     *
     * \param[in] global Coordinates in the reference triangle
     * \returns Number of the subtriangle containing <tt>global</tt>
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
       \param[out] subElement Number of the subtriangle containing <tt>global</tt>
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
     \brief Base class for LocalBasis classes based on uniform refinement in 3D; provides numbering and local coordinates of subelements

     Shape functions like these are necessary for hierarchical error estimators
     for certain nonlinear problems.

     \tparam D Type to represent the field in the domain.

     \nosubgrouping
   */
  template<class D>
  class RefinedSimplexLocalBasis<D,3>
  {
  protected:

    /** \brief Protected default constructor so this class can only be instantiated as a base class. */
    RefinedSimplexLocalBasis() {}

    /** \brief Get the number of the subsimplex containing a given point in the reference element
     *
     * Defining the following points in the reference simplex
     *
     * 0: (0.0, 0.0, 0.0)
     * 1: (1.0, 0.0, 0.0)
     * 2: (0.0, 1.0, 0.0)
     * 3: (0.0, 0.0, 1.0)
     * 4: (0.5, 0.0, 0.0)
     * 5: (0.5, 0.5, 0.0)
     * 6: (0.0, 0.5, 0.0)
     * 7: (0.0, 0.0, 0.5)
     * 8: (0.5, 0.0, 0.5)
     * 9: (0.0, 0.5, 0.5)
     *
     * The subsimplices are numbered according to
     *
     * 0: 0467  -
     * 1: 4158   |_ "cut off" vertices
     * 2: 6529   |
     * 3: 7893  -
     *
     * 4: 6487  -
     * 5: 4568   |_  octahedron partition
     * 6: 6897   |
     * 7: 6895  -
     *
     * \param[in] global Coordinates in the reference simplex
     * \returns Number of the subsimplex containing <tt>global</tt>
     */
    static int getSubElement(const FieldVector<D,3>& global)
    {
      if (global[0] + global[1] + global[2] <= 0.5)
        return 0;
      else if (global[0] >= 0.5)
        return 1;
      else if (global[1] >= 0.5)
        return 2;
      else if (global[2] >= 0.5)
        return 3;
      else if ((global[0] + global[1] <= 0.5)and (global[1] + global[2] <= 0.5))
        return 4;
      else if ((global[0] + global[1] >= 0.5)and (global[1] + global[2] <= 0.5))
        return 5;
      else if ((global[0] + global[1] <= 0.5)and (global[1] + global[2] >= 0.5))
        return 6;
      else if ((global[0] + global[1] >= 0.5)and (global[1] + global[2] >= 0.5))
        return 7;

      DUNE_THROW(InvalidStateException, "no subelement defined");

    }
    /** \brief Get local coordinates in the subsimplex

       \param[in] global Coordinates in the reference simplex
       \param[out] subElement Number of the subsimplex containing <tt>global</tt>
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
