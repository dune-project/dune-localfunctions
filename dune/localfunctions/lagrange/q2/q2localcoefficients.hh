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

        li[ 0] = LocalKey(0,3,0);
        li[ 1] = LocalKey(6,2,0);
        li[ 2] = LocalKey(1,3,0);
        li[ 3] = LocalKey(4,2,0);
        li[ 4] = LocalKey(4,1,0);
        li[ 5] = LocalKey(5,2,0);
        li[ 6] = LocalKey(2,3,0);
        li[ 7] = LocalKey(7,2,0);
        li[ 8] = LocalKey(3,3,0);

        li[ 9] = LocalKey(0,2,0);
        li[10] = LocalKey(2,1,0);
        li[11] = LocalKey(1,2,0);
        li[12] = LocalKey(0,1,0);
        li[13] = LocalKey(0,0,0);
        li[14] = LocalKey(1,1,0);
        li[15] = LocalKey(2,2,0);
        li[16] = LocalKey(3,1,0);
        li[17] = LocalKey(3,2,0);

        li[18] = LocalKey(4,3,0);
        li[19] = LocalKey(10,2,0);
        li[20] = LocalKey(5,3,0);
        li[21] = LocalKey(8,2,0);
        li[22] = LocalKey(5,1,0);
        li[23] = LocalKey(9,2,0);
        li[24] = LocalKey(6,3,0);
        li[25] = LocalKey(11,2,0);
        li[26] = LocalKey(7,3,0);

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
