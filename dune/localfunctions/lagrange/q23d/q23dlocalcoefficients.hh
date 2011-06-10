// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_Q23DLOCALCOEFFICIENTS_HH
#define DUNE_Q23DLOCALCOEFFICIENTS_HH

#include <cstddef>
#include <iostream>
#include <vector>

#include <dune/localfunctions/common/localkey.hh>

namespace Dune
{

  /**@ingroup LocalLayoutImplementation
           \brief Layout map for Q2 elements

           \nosubgrouping
       \implements Dune::LocalCoefficientsVirtualImp
   */
  class Q23DLocalCoefficients
  {
  public:
    //! \brief Standard constructor
    Q23DLocalCoefficients () : li(27)
    {
      for (std::size_t i = 0; i < 8; i++)
        li[i] = LocalKey(i, 3, 0);

      for (std::size_t i = 8; i < 20; i++)
        li[i] = LocalKey(i-8, 2, 0);

      for (std::size_t i = 20; i < 26; i++)
        li[i] = LocalKey(i-20, 1, 0);

      li[26] = LocalKey(0, 0, 0);
    }

    //! number of coefficients
    std::size_t size () const
    {
      return 27;
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
