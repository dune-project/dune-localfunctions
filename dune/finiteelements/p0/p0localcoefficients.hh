// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_P0LOCALCOEFFICIENTS_HH
#define DUNE_P0LOCALCOEFFICIENTS_HH

#include <iostream>
#include <vector>

#include "../common/localcoefficients.hh"

namespace Dune
{

  /**@ingroup LocalLayoutImplementation
         \brief Layout map for P0 elements

         \nosubgrouping
   */
  class P0LocalCoefficients : public LocalCoefficientsInterface<P0LocalCoefficients>
  {
  public:
    //! \brief Standard constructor
    P0LocalCoefficients () : index(0,0,0)
    {}

    //! number of coefficients
    int size ()
    {
      return 1;
    }

    //! get i'th index
    const LocalIndex& localIndex (int i) const
    {
      return index;
    }

  private:
    LocalIndex index;
  };

}
#endif
