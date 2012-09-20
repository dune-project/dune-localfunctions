// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_LOCALFUNCTIONS_QK_LOCALCOEFFICIENTS_HH
#define DUNE_LOCALFUNCTIONS_QK_LOCALCOEFFICIENTS_HH


namespace Dune
{
  template<int k, int d>
  class QkLocalCoefficients {       };


  template<int d>
  class QkLocalCoefficients<1,d>
  {
    enum { k = 1 };
  public:
    //! \brief Standard constructor
    QkLocalCoefficients () : li(Power_m_p<k+1,d>::power)
    {
      for (std::size_t i=0; i<Power_m_p<k+1,d>::power; i++)
        li[i] = LocalKey(i,d,0);
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

  template<>
  class QkLocalCoefficients<2,2>
  {
    enum { k = 2 };
    enum { d = 2 };
  public:
    //! \brief Standard constructor
    QkLocalCoefficients () : li(Power_m_p<k+1,d>::power)
    {
      li[0] = LocalKey(0,2,0);
      li[1] = LocalKey(2,1,0);
      li[2] = LocalKey(1,2,0);
      li[3] = LocalKey(0,1,0);
      li[4] = LocalKey(0,0,0);
      li[5] = LocalKey(1,1,0);
      li[6] = LocalKey(2,2,0);
      li[7] = LocalKey(3,1,0);
      li[8] = LocalKey(3,2,0);
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

  /** Unfortunately this specialization is needed to avoid that the compiler has to decide between <1,d> and <k,2> (which he can't) **/
  template<>
  class QkLocalCoefficients<1,2>
  {
    enum { k = 1 };
    enum { d = 2 };
  public:
    //! \brief Standard constructor
    QkLocalCoefficients () : li(Power_m_p<k+1,d>::power)
    {
      for (std::size_t i=0; i<Power_m_p<k+1,d>::power; i++)
        li[i] = LocalKey(i,d,0);
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

  /** \todo Please doc me !
   * arbritary order in two dimensions
   */
  template<int k>
  class QkLocalCoefficients<k,2>
  {

    enum { d = 2 };
  public:
    //! \brief Standard constructor
    QkLocalCoefficients () : li(Power_m_p<k+1,d>::power)
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
      li[lastIndex++] = LocalKey(0,2,0);    // corner 0
      for (unsigned i = 0; i < k - 1; ++i) {
        li[lastIndex++] = LocalKey(2,1,i);   // inner dofs of lower edge (2)
      }
      li[lastIndex++] = LocalKey(1,2,0);   // corner 1

      // iterate from bottom to top over inner edge dofs
      for (unsigned e = 0; e < k - 1; ++e) {
        li[lastIndex++] = LocalKey(0,1,e);   // left edge (0)
        for (unsigned i = 0; i < k - 1; ++i) {
          li[lastIndex++] = LocalKey(0,0,lastInnerFaceIndex++);   // face dofs
        }
        li[lastIndex++] = LocalKey(1,1,e);   // right edge (1)
      }

      // upper edge (3)
      li[lastIndex++] = LocalKey(2,2,0);   // corner 2
      for (unsigned i = 0; i < k - 1; ++i) {
        li[lastIndex++] = LocalKey(3,1,i);   // inner dofs of upper edge (3)
      }
      li[lastIndex++] = LocalKey(3,2,0);   // corner 3

      const unsigned numIndices = Power_m_p<k+1,d>::power;

      --lastIndex; --lastInnerFaceIndex;

      assert(numIndices==lastIndex+1);
      assert(((k-1)*(k-1))==lastInnerFaceIndex+1);
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

  //// d = 3

  /** Unfortunately this specialization is needed to avoid that the compiler has to decide between <1,d> and <k,2> (which he can't) **/
  template<>
  class QkLocalCoefficients<1,3>
  {
    enum { k = 1 };
    enum { d = 3 };
  public:
    //! \brief Standard constructor
    QkLocalCoefficients () : li(Power_m_p<k+1,d>::power)
    {
      for (std::size_t i=0; i<Power_m_p<k+1,d>::power; i++)
        li[i] = LocalKey(i,d,0);
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

  template<int k>
  class QkLocalCoefficients<k,3>
  {
    enum { d = 3 };
  public:
    //! \brief Standard constructor
    QkLocalCoefficients () : li(Power_m_p<k+1,d>::power)
    {
      assert(k>0);
      unsigned lastIndex=0;
      unsigned lastInnerFaceIndex=0;
      const unsigned numIndices = Power_m_p<k+1,d>::power;
      const unsigned numFaceIndices = Power_m_p<k+1,d-1>::power;
      const unsigned numInnerEdgeDofs = k-1;
      const unsigned numInnerFaceDofs = numInnerEdgeDofs * numInnerEdgeDofs;
      // const unsigned numInnerVolumeDofs = numInnerEdgeDofs * numInnerEdgeDofs *numInnerEdgeDofs;

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
      li[lastIndex++] = LocalKey(0,d,0);   // corner 0
      for (unsigned i = 0; i < numInnerEdgeDofs; ++i) {
        li[lastIndex++] = LocalKey(6,d-1,i);   // inner dofs of lower edge (2)
      }
      li[lastIndex++] = LocalKey(1,d,0);   // corner 1

      // iterate from bottom to top over inner edge dofs
      for (unsigned e = 0; e < numInnerEdgeDofs; ++e) {
        li[lastIndex++] = LocalKey(4,d-1,e);   // left edge (4)
        for (unsigned i = 0; i < numInnerEdgeDofs; ++i) {
          li[lastIndex++] = LocalKey(4,d-2,lastInnerFaceIndex++);   // inner face dofs
        }
        li[lastIndex++] = LocalKey(5,d-1,e);   // right edge (5)
      }

      // upper edge (7)
      li[lastIndex++] = LocalKey(2,d,0);   // corner 2
      for (unsigned i = 0; i < k - 1; ++i) {
        li[lastIndex++] = LocalKey(7,d-1,i);   // inner dofs of upper edge (7)
      }
      li[lastIndex++] = LocalKey(3,d,0);   // corner 3

      assert(((k-1)*(k-1))==lastInnerFaceIndex);

      assert(numFaceIndices==lastIndex);   // added 1 face so far
      /////////////////////////////////////////// end bottom face (4)

      ///////////////////// inner faces
      for(unsigned f = 0; f < numInnerEdgeDofs; ++f) {

        lastInnerFaceIndex=0;

        // lower edge (connecting  edges 0 and 1)
        li[lastIndex++] = LocalKey(0,d-1,f);   // dof on edge 0
        for (unsigned i = 0; i < numInnerEdgeDofs; ++i) {
          li[lastIndex++] = LocalKey(2,d-2,i+(f*numInnerEdgeDofs));   // dof in front face (2)
        }
        li[lastIndex++] = LocalKey(1,d-1,f);   // dof on edge 1

        // iterate from bottom to top over inner edge dofs
        for (unsigned e = 0; e < numInnerEdgeDofs; ++e) {
          li[lastIndex++] = LocalKey(0,d-2,e+f*numInnerEdgeDofs);   // on left face (0)
          for (unsigned i = 0; i < numInnerEdgeDofs; ++i) {
            li[lastIndex++] = LocalKey(0,d-3,i + f*numInnerFaceDofs + e*numInnerEdgeDofs);   // volume dofs
          }
          li[lastIndex++] = LocalKey(1,d-2,e+f*numInnerEdgeDofs);   // right edge (5)
        }

        // upper edge (connecting  edges 0 and 1)
        li[lastIndex++] = LocalKey(2,d-1,f);   // dof on edge 2
        for (unsigned i = 0; i < numInnerEdgeDofs; ++i) {
          li[lastIndex++] = LocalKey(3,d-2,i+(f*numInnerEdgeDofs));   // dof in rear face (3)
        }
        li[lastIndex++] = LocalKey(3,d-1,f);   // dof on edge 3

        assert(lastIndex==(f+1+1)*numFaceIndices);
      }

      ////////////////////////////////////////// top face (5)
      lastInnerFaceIndex=0;
      // lower edge (10)
      li[lastIndex++] = LocalKey(4,d,0);   // corner 4
      for (unsigned i = 0; i < k - 1; ++i) {
        li[lastIndex++] = LocalKey(10,d-1,i);   // inner dofs of lower edge (10)
      }
      li[lastIndex++] = LocalKey(5,d,0);   // corner 5

      // iterate from bottom to top over inner edge dofs
      for (unsigned e = 0; e < k - 1; ++e) {
        li[lastIndex++] = LocalKey(8,d-1,e);   // left edge (8)
        for (unsigned i = 0; i < k - 1; ++i) {
          li[lastIndex++] = LocalKey(5,d-2,lastInnerFaceIndex++);   // face dofs
        }
        li[lastIndex++] = LocalKey(9,d-1,e);   // right edge (9)
      }

      // upper edge (11)
      li[lastIndex++] = LocalKey(6,d,0);   // corner 6
      for (unsigned i = 0; i < k - 1; ++i) {
        li[lastIndex++] = LocalKey(11,d-1,i);   // inner dofs of upper edge (7)
      }
      li[lastIndex++] = LocalKey(7,d,0);   // corner 7

      assert(numIndices==lastIndex);
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
