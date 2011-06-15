// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_Q22DLOCALCOEFFICIENTS_HH
#define DUNE_Q22DLOCALCOEFFICIENTS_HH

#warning This file is deprecated and will be removed after Dune 2.2.\
  Please use q2localcoefficients.hh instead!

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
  class Q22DLocalCoefficients
  {
  public:
    //! \brief Standard constructor
    Q22DLocalCoefficients () : li(9)
    {
      for (std::size_t i=0; i<9; i++)
        li[i] = LocalKey(i%4,2-i/4,0);
    }

    //! number of coefficients
    std::size_t size () const
    {
      return 9;
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
