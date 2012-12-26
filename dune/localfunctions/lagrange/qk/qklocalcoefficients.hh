// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_LOCALFUNCTIONS_QK_LOCALCOEFFICIENTS_HH
#define DUNE_LOCALFUNCTIONS_QK_LOCALCOEFFICIENTS_HH

#include <cassert>
#include <vector>

#include <dune/common/exceptions.hh>
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

    // Return i as a d-digit number in the (k+1)-nary system
    static array<unsigned int,d> multiindex (unsigned int i)
    {
      array<unsigned int,d> alpha;
      for (int j=0; j<d; j++)
      {
        alpha[j] = i % (k+1);
        i = i/(k+1);
      }
      return alpha;
    }


    void setup2d(std::vector<unsigned int>& subEntity,
                 std::vector<unsigned int>& index)
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
      subEntity[lastIndex] = 0;                 // corner 0
      index[lastIndex++] = 0;
      for (unsigned i = 0; i < k - 1; ++i) {
        subEntity[lastIndex] = 2;           // inner dofs of lower edge (2)
        index[lastIndex++] = i;
      }

      subEntity[lastIndex] = 1;                 // corner 1
      index[lastIndex++] = 0;

      // iterate from bottom to top over inner edge dofs
      for (unsigned e = 0; e < k - 1; ++e) {
        subEntity[lastIndex] = 0;                   // left edge (0)
        index[lastIndex++] = e;
        for (unsigned i = 0; i < k - 1; ++i) {
          subEntity[lastIndex] = 0;                     // face dofs
          index[lastIndex++] = lastInnerFaceIndex++;
        }
        subEntity[lastIndex] = 1;                   // right edge (1)
        index[lastIndex++] = e;
      }

      // upper edge (3)
      subEntity[lastIndex] = 2;                 // corner 2
      index[lastIndex++] = 0;
      for (unsigned i = 0; i < k - 1; ++i) {
        subEntity[lastIndex] = 3;                   // inner dofs of upper edge (3)
        index[lastIndex++] = i;
      }

      subEntity[lastIndex] = 3;                 // corner 3
      index[lastIndex++] = 0;

      assert((StaticPower<k+1,d>::power==lastIndex));
      assert(((k-1)*(k-1))==lastInnerFaceIndex);
    }


    void setup3d()
    {
      assert(k>0);
      unsigned lastIndex=0;
      unsigned lastInnerFaceIndex=0;
#ifndef NDEBUG
      const unsigned numIndices = StaticPower<k+1,d>::power;
      const unsigned numFaceIndices = StaticPower<k+1,d-1>::power;
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
    QkLocalCoefficients () : li(StaticPower<k+1,d>::power)
    {
      // Set up array of codimension-per-dof-number
      std::vector<unsigned int> codim(li.size());

      for (std::size_t i=0; i<codim.size(); i++) {
        codim[i] = 0;
        // Codimension gets reduced by 1 for each coordinate direction
        // where dof is on boundary
        array<unsigned int,d> mIdx = multiindex(i);
        for (int j=0; j<d; j++)
          if (mIdx[j]==0 or mIdx[j]==k)
            codim[i]--;
      }

      // Set up entity and dof numbers for each (supported) dimension separately
      std::vector<unsigned int> subEntity(li.size());
      std::vector<unsigned int> index(li.size());

      if (k==1) {

        for (std::size_t i=0; i<StaticPower<k+1,d>::power; i++)
          li[i] = LocalKey(i,codim[i],0);

      } else if (d==2) {

        setup2d(subEntity, index);

        for (size_t i=0; i<li.size(); i++)
          li[i] = LocalKey(subEntity[i], codim[i], index[i]);

      } else if (d==3) {

        setup3d();

      } else
        DUNE_THROW(Dune::NotImplemented, "QkLocalCoefficients for k==" << k << " and d==" << d);
    }

    //! number of coefficients
    std::size_t size () const
    {
      return StaticPower<k+1,d>::power;
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
