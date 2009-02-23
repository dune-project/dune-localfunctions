// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_Q13DLOCALCOEFFICIENTS_HH
#define DUNE_Q13DLOCALCOEFFICIENTS_HH

#include <iostream>
#include <vector>

#include "../common/localcoefficients.hh"

namespace Dune
{

  /**@ingroup LocalLayoutImplementation
         \brief Layout map for Q1 elements

         \nosubgrouping
   */
  class Q13DLocalCoefficients : public LocalCoefficientsInterface<Q13DLocalCoefficients>
  {
  public:
    //! \brief Standard constructor
    Q13DLocalCoefficients () : li(8)
    {
      for (int i=0; i<8; i++)
        li[i] = LocalIndex(i,3,0);
    }

    //! number of coefficients
    int size () const
    {
      return 8;
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
