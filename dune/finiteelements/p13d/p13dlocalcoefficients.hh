// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_P13DLOCALCOEFFICIENTS_HH
#define DUNE_P13DLOCALCOEFFICIENTS_HH

#include <iostream>
#include <vector>

#include "../common/localcoefficients.hh"

namespace Dune
{

  /**@ingroup LocalLayoutImplementation
         \brief Layout map for P13D elements

         \nosubgrouping
   */
  class P13DLocalCoefficients : public LocalCoefficientsInterface<P13DLocalCoefficients>
  {
  public:
    //! \brief Standard constructor
    P13DLocalCoefficients () : li(4)
    {
      for (int i=0; i<4; i++)
        li[i] = LocalIndex(i,2,0);
    }

    //! number of coefficients
    int size () const
    {
      return 4;
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
