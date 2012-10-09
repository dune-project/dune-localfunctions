// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_LOCALFUNCTIONS_QK_LOCALCOEFFICIENTS_HH
#define DUNE_LOCALFUNCTIONS_QK_LOCALCOEFFICIENTS_HH

#include <cassert>

#include <dune/common/power.hh>

#include <dune/localfunctions/common/localkey.hh>

namespace Dune
{
  /** \brief Attaches a shape function to an entity
   *
   * \tparam k Polynomial order
   * \tparam d Dimension of the reference cube
   */
  template<int k, int d>
  class QkLocalCoefficients {

    void setup2d()
    {
      assert(k>0);
      unsigned lastIndex=0;
      unsigned lastInnerFaceIndex=0;

      // LocalKey: entity number , entity codim, dof indices within each entity
      /* edge and vertex numbering
         2----3----3
       |         |
       |         |
         0         1
       |         |
       |         |
         0----2----1
       */

      // lower edge (2)
      li[lastIndex++] = LocalKey(0,2,0);        // corner 0
      for (unsigned i = 0; i < k - 1; ++i) {
        li[lastIndex++] = LocalKey(2,1,i);         // inner dofs of lower edge (2)
      }
      li[lastIndex++] = LocalKey(1,2,0);       // corner 1

      // iterate from bottom to top over inner edge dofs
      for (unsigned e = 0; e < k - 1; ++e) {
        li[lastIndex++] = LocalKey(0,1,e);         // left edge (0)
        for (unsigned i = 0; i < k - 1; ++i)
          li[lastIndex++] = LocalKey(0,0,lastInnerFaceIndex++);           // face dofs

        li[lastIndex++] = LocalKey(1,1,e);         // right edge (1)
      }

      // upper edge (3)
      li[lastIndex++] = LocalKey(2,2,0);       // corner 2
      for (unsigned i = 0; i < k - 1; ++i) {
        li[lastIndex++] = LocalKey(3,1,i);         // inner dofs of upper edge (3)
      }
      li[lastIndex++] = LocalKey(3,2,0);       // corner 3

#ifndef NDEBUG
      const unsigned numIndices = Power_m_p<k+1,d>::power;
#endif

      --lastIndex; --lastInnerFaceIndex;

      assert(numIndices==lastIndex+1);
      assert(((k-1)*(k-1))==lastInnerFaceIndex+1);
    }


    void setup3d()
    {
      assert(k>0);
      unsigned lastIndex=0;
      unsigned lastInnerFaceIndex=0;
#ifndef NDEBUG
      const unsigned numIndices = Power_m_p<k+1,d>::power;
      const unsigned numFaceIndices = Power_m_p<k+1,d-1>::power;
#endif
      const unsigned numInnerEdgeDofs = k-1;
      const unsigned numInnerFaceDofs = numInnerEdgeDofs * numInnerEdgeDofs;

      // LocalKey: entity number , entity codim, dof indices within each entity
      /* edge and vertex numbering
          2----3----3
       |         |
       |         |
          0         1
       |         |
       |         |
          0----2----1
       */

      // bottom face (4)
      lastIndex=0;
      lastInnerFaceIndex=0;
      // lower edge (2)
      li[lastIndex++] = LocalKey(0,d,0);       // corner 0
      for (unsigned i = 0; i < numInnerEdgeDofs; ++i)
        li[lastIndex++] = LocalKey(6,d-1,i);         // inner dofs of lower edge (2)

      li[lastIndex++] = LocalKey(1,d,0);       // corner 1

      // iterate from bottom to top over inner edge dofs
      for (unsigned e = 0; e < numInnerEdgeDofs; ++e) {
        li[lastIndex++] = LocalKey(4,d-1,e);         // left edge (4)
        for (unsigned i = 0; i < numInnerEdgeDofs; ++i) {
          li[lastIndex++] = LocalKey(4,d-2,lastInnerFaceIndex++);           // inner face dofs
        }
        li[lastIndex++] = LocalKey(5,d-1,e);         // right edge (5)
      }

      // upper edge (7)
      li[lastIndex++] = LocalKey(2,d,0);       // corner 2
      for (unsigned i = 0; i < k - 1; ++i) {
        li[lastIndex++] = LocalKey(7,d-1,i);         // inner dofs of upper edge (7)
      }
      li[lastIndex++] = LocalKey(3,d,0);       // corner 3

      assert(((k-1)*(k-1))==lastInnerFaceIndex);

      assert(numFaceIndices==lastIndex);       // added 1 face so far
      /////////////////////////////////////////// end bottom face (4)

      ///////////////////// inner faces
      for(unsigned f = 0; f < numInnerEdgeDofs; ++f) {

        lastInnerFaceIndex=0;

        // lower edge (connecting  edges 0 and 1)
        li[lastIndex++] = LocalKey(0,d-1,f);         // dof on edge 0
        for (unsigned i = 0; i < numInnerEdgeDofs; ++i) {
          li[lastIndex++] = LocalKey(2,d-2,i+(f*numInnerEdgeDofs));           // dof in front face (2)
        }
        li[lastIndex++] = LocalKey(1,d-1,f);         // dof on edge 1

        // iterate from bottom to top over inner edge dofs
        for (unsigned e = 0; e < numInnerEdgeDofs; ++e) {
          li[lastIndex++] = LocalKey(0,d-2,e+f*numInnerEdgeDofs);           // on left face (0)
          for (unsigned i = 0; i < numInnerEdgeDofs; ++i) {
            li[lastIndex++] = LocalKey(0,d-3,i + f*numInnerFaceDofs + e*numInnerEdgeDofs);             // volume dofs
          }
          li[lastIndex++] = LocalKey(1,d-2,e+f*numInnerEdgeDofs);           // right edge (5)
        }

        // upper edge (connecting  edges 0 and 1)
        li[lastIndex++] = LocalKey(2,d-1,f);         // dof on edge 2
        for (unsigned i = 0; i < numInnerEdgeDofs; ++i) {
          li[lastIndex++] = LocalKey(3,d-2,i+(f*numInnerEdgeDofs));           // dof in rear face (3)
        }
        li[lastIndex++] = LocalKey(3,d-1,f);         // dof on edge 3

        assert(lastIndex==(f+1+1)*numFaceIndices);
      }

      ////////////////////////////////////////// top face (5)
      lastInnerFaceIndex=0;
      // lower edge (10)
      li[lastIndex++] = LocalKey(4,d,0);       // corner 4
      for (unsigned i = 0; i < k - 1; ++i) {
        li[lastIndex++] = LocalKey(10,d-1,i);         // inner dofs of lower edge (10)
      }
      li[lastIndex++] = LocalKey(5,d,0);       // corner 5

      // iterate from bottom to top over inner edge dofs
      for (unsigned e = 0; e < k - 1; ++e) {
        li[lastIndex++] = LocalKey(8,d-1,e);         // left edge (8)
        for (unsigned i = 0; i < k - 1; ++i) {
          li[lastIndex++] = LocalKey(5,d-2,lastInnerFaceIndex++);           // face dofs
        }
        li[lastIndex++] = LocalKey(9,d-1,e);         // right edge (9)
      }

      // upper edge (11)
      li[lastIndex++] = LocalKey(6,d,0);       // corner 6
      for (unsigned i = 0; i < k - 1; ++i)
        li[lastIndex++] = LocalKey(11,d-1,i);             // inner dofs of upper edge (7)

      li[lastIndex++] = LocalKey(7,d,0);       // corner 7

      assert(numIndices==lastIndex);
    }

  public:
    //! \brief Default constructor
    QkLocalCoefficients () : li(Power_m_p<k+1,d>::power)
    {
      if (k==1) {

        for (std::size_t i=0; i<Power_m_p<k+1,d>::power; i++)
          li[i] = LocalKey(i,d,0);

      } else if (d==2) {

        setup2d();

      } else if (d==3) {

        setup3d();

      } else
        DUNE_THROW(Dune::NotImplemented, "QkLocalCoefficients for k==" << k << " and d==" << d);
    }

    //! number of coefficients
    std::size_t size () const
    {
      return Power_m_p<k+1,d>::power;
    }

    //! get i'th index
    const LocalKey& localKey (std::size_t i) const
    {
      return li[i];
    }

  private:
    std::vector<LocalKey> li;
  };

}

#endif
