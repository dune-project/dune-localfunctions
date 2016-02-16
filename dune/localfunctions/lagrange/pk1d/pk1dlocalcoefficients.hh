// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PK1DLOCALCOEFFICIENTS_HH
#define DUNE_PK1DLOCALCOEFFICIENTS_HH

#include <cstddef>
#include <vector>

#include <dune/localfunctions/common/localkey.hh>

namespace Dune
{

  /**@ingroup LocalLayoutImplementation
         \brief Layout map for Pk elements

         \nosubgrouping
     \implements Dune::LocalCoefficientsVirtualImp
   */
  template<unsigned int k>
  class Pk1DLocalCoefficients
  {
    enum {N = k+1};

  public:
    //! \brief Standard constructor
    Pk1DLocalCoefficients () : li(N)
    {
      fill_default();
    }

    //! constructor for eight variants with order on edges flipped
    Pk1DLocalCoefficients (int variant) : li(N)
    {
      fill_default();
    }

    /** Constructor for six variants with permuted vertices.

        \param vertexmap The permutation of the vertices.  This
        can for instance be generated from the global indices of
        the vertices by reducing those to the integers 0...2.  This may be any
        object which for which the expression \c vertexmap[i] is defined
        appropriately (like an array, a pointer, a std::vector, or a
        random-access iterator.
     */
    template<class VertexMap>
    explicit Pk1DLocalCoefficients(const VertexMap &vertexmap) : li(N)
    {
      fill_default();
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

    void fill_default()
    {
      li.resize(N);

      if (N==1) {
        li[0] = LocalKey(0,0,0);
      } else {
        li[0] = LocalKey(0,1,0);
        for (int i=1; i<N-1; i++)
          li[i] = LocalKey(0,0,i-1);            // element dofs
        li.back() = LocalKey(1,1,0);
      }
    }
  };

}

#endif
