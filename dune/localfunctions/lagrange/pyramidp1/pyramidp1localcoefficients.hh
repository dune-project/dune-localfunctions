// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PYRAMID_P1_LOCALCOEFFICIENTS_HH
#define DUNE_PYRAMID_P1_LOCALCOEFFICIENTS_HH

#include <cstddef>
#include <iostream>
#include <vector>

#include <dune/localfunctions/common/localkey.hh>

namespace Dune
{

  /**@ingroup LocalLayoutImplementation
         \brief Layout map for PyramidP1 elements

         \nosubgrouping
   */
  class PyramidP1LocalCoefficients
  {
  public:
    //! \brief Standard constructor
    PyramidP1LocalCoefficients () : li(5)
    {
      for (std::size_t i=0; i<5; i++)
        li[i] = LocalKey(i,3,0);
    }

    //! number of coefficients
    std::size_t size () const
    {
      return 5;
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
