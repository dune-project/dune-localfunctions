// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_EDGES12DLOCALCOEFFICIENTS_HH
#define DUNE_EDGES12DLOCALCOEFFICIENTS_HH

#include <iostream>
#include <vector>

#include "../common/localcoefficients.hh"

namespace Dune
{

  /**@ingroup LocalLayoutImplementation
     \brief Layout map for P0 elements

     \nosubgrouping
   */
  class EdgeS12DLocalCoefficients
#ifdef DUNE_VIRTUAL_SHAPEFUNCTIONS
    : public LocalCoefficientsInterface
#else
    : public LocalCoefficientsInterface<EdgeS12DLocalCoefficients>
#endif
  {
  public:
    //! \brief Standard constructor
    EdgeS12DLocalCoefficients () : li(3)
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

#endif // DUNE_EDGES12DLOCALCOEFFICIENTS_HH
