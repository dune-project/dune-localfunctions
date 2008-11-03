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
  class RT0TriangleLocalCoefficients : public LocalCoefficientsInterface<RT0TriangleLocalCoefficients>
  {
  public:
    //! \brief Standard constructor
    RT0TriangleLocalCoefficients () : li(3)
    {
      for (int i=0; i<3; i++)
        li[i] = LocalIndex(i,1,0);
    }

    //! number of coefficients
    int size () const
    {
      return 3;
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
