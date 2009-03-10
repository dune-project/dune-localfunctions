// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_RT0TRIANGLELOCALCOEFFICIENTS_HH
#define DUNE_RT0TRIANGLELOCALCOEFFICIENTS_HH

#include <iostream>
#include <vector>

#include "../common/localcoefficients.hh"

namespace Dune
{

  /**@ingroup LocalLayoutImplementation
         \brief Layout map for RT0 elements

         \nosubgrouping
   */
  class RT02DLocalCoefficients : public LocalCoefficientsInterface<RT02DLocalCoefficients>
  {
  public:
    //! \brief Standard constructor
    RT02DLocalCoefficients () : li(3)
    {
      for (int i=0; i<3; i++)
        li[i] = LocalKey(i,1,0);
    }

    //! number of coefficients
    int size () const
    {
      return 3;
    }

    //! get i'th index
    const LocalKey& localKey (int i) const
    {
      return li[i];
    }

  private:
    std::vector<LocalKey> li;
  };

}

#endif
