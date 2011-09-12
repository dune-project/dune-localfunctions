// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_Q2_LOCALCOEFFICIENTS_HH
#define DUNE_Q2_LOCALCOEFFICIENTS_HH

#include <cstddef>
#include <vector>

#include <dune/localfunctions/common/localkey.hh>

namespace Dune
{

  /**@ingroup LocalLayoutImplementation
           \brief Layout map for Q2 elements

      \tparam dim Dimension of the reference element

           \nosubgrouping
       \implements Dune::LocalCoefficientsVirtualImp
   */
  template <int dim>
  class Q2LocalCoefficients
  {
  public:
    //! \brief Standard constructor
    Q2LocalCoefficients () : li(size())
    {
      switch (dim) {

      case 1 : {
        li[0] = LocalKey(0,1,0);       // left vertex
        li[1] = LocalKey(0,0,0);       // element
        li[2] = LocalKey(1,1,0);       // right vertex
        break;
      }

      case 2 : {

        li[0] = LocalKey(0,2,0);
        li[1] = LocalKey(2,1,0);
        li[2] = LocalKey(1,2,0);
        li[3] = LocalKey(0,1,0);
        li[4] = LocalKey(0,0,0);
        li[5] = LocalKey(1,1,0);
        li[6] = LocalKey(2,2,0);
        li[7] = LocalKey(3,1,0);
        li[8] = LocalKey(3,2,0);

        break;
      }

      case 3 : {

        for (std::size_t i = 0; i < 8; i++)
          li[i] = LocalKey(i, 3, 0);

        for (std::size_t i = 8; i < 20; i++)
          li[i] = LocalKey(i-8, 2, 0);

        for (std::size_t i = 20; i < 26; i++)
          li[i] = LocalKey(i-20, 1, 0);

        li[26] = LocalKey(0, 0, 0);

        break;
      }
      default :
        DUNE_THROW(NotImplemented, "Q2LocalCoefficients for dim==" << dim);
      }
    }

    //! number of coefficients
    std::size_t size () const
    {
      int size = 1;
      for (int i=0; i<dim; i++)
        size *= 3;
      return size;
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
