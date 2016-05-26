// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PK3DLOCALCOEFFICIENTS_HH
#define DUNE_PK3DLOCALCOEFFICIENTS_HH

#include <cstddef>
#include <iostream>
#include <vector>

#include <dune/localfunctions/common/localkey.hh>

namespace Dune
{

  /**@ingroup LocalLayoutImplementation
     \brief Please doc me!

     \nosubgrouping
     \implements Dune::LocalCoefficientsVirtualImp
   */
  template<unsigned int k>
  class Pk3DLocalCoefficients
  {
    enum {N = (k+1)*(k+2)*(k+3)/6};

  public:
    //! \brief Standard constructor
    Pk3DLocalCoefficients () : li(N)
    {
      const unsigned int vertexmap[4] = {0, 1, 2, 3};
      generate_local_keys(vertexmap);
    }

    /** Constructor for variants with permuted vertices.

        \param vertexmap The permutation of the vertices.  This
        can for instance be generated from the global indices of
        the vertices by reducing those to the integers 0...3
     */
    Pk3DLocalCoefficients (const unsigned int vertexmap[4]) : li(N)
    {
      generate_local_keys(vertexmap);
    }

    //! number of coefficients
    std::size_t size () const
    {
      return N;
    }

    //! get i'th index
    const LocalKey& localKey (std::size_t i) const
    {
      return li[i];
    }

  private:
    std::vector<LocalKey> li;

    void generate_local_keys(const unsigned int vertexmap[4])
    {
      unsigned int subindex[16];
      unsigned int codim_count[4] = {0};
      for (unsigned int m = 1; m < 16; ++m)
      {
        unsigned int codim = !(m&1) + !(m&2) + !(m&4) + !(m&8);
        subindex[m] = codim_count[codim]++;
      }

      int a1 = (3*k + 12)*k + 11;
      int a2 = -3*k - 6;
      unsigned int dof_count[16] = {0};
      unsigned int i[4];
      for (i[3] = 0; i[3] <= k; ++i[3])
        for (i[2] = 0; i[2] <= k - i[3]; ++i[2])
          for (i[1] = 0; i[1] <= k - i[2] - i[3]; ++i[1])
          {
            i[0] = k - i[1] - i[2] - i[3];
            unsigned int j[4];
            unsigned int entity = 0;
            unsigned int codim = 0;
            for (unsigned int m = 0; m < 4; ++m)
            {
              j[m] = i[vertexmap[m]];
              entity += !!j[m] << m;
              codim += !j[m];
            }
            int local_index = j[3]*(a1 + (a2 + j[3])*j[3])/6
                              + j[2]*(2*(k - j[3]) + 3 - j[2])/2 + j[1];
            li[local_index] = LocalKey(subindex[entity], codim, dof_count[entity]++);
          }
    }
  };

}

#endif
