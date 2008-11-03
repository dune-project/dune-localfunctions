// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_Q22DLOCALCOEFFICIENTS_HH
#define DUNE_Q22DLOCALCOEFFICIENTS_HH

#include <iostream>
#include <vector>

#include "locallayout.hh"

#include "../common/localcoefficients.hh"

namespace Dune
{

  /**@ingroup LocalLayoutImplementation
         \brief Layout map for Q2 elements

         \nosubgrouping
   */
  class Q22DLocalCoefficients : public LocalCoefficientsInterface<Q22DLocalCoefficients>
  {
  public:
    //! \brief Standard constructor
    Q22DLocalCoefficients () : li(9)
    {
      for (int i=0; i<9; i++)
        li[i] = LocalIndex(i%4,2-i/4,0);
    }

    //! number of coefficients
    int size ()
    {
      return 9;
    }

    //! get i'th index
    const LocalIndex& localIndex (int i) const
    {
      return li[i];
    }

  private:
    std::vector<LocalIndex> li;
  };

}

#endif
