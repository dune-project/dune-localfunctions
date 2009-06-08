// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_EDGES03DLOCALCOEFFICIENTS_HH
#define DUNE_EDGES03DLOCALCOEFFICIENTS_HH

#include <iostream>
#include <vector>

#include "../common/localcoefficients.hh"

namespace Dune
{

  /**@ingroup LocalLayoutImplementation
     \brief Layout map for lowest order edge elements on tetrahedrons

     \nosubgrouping
   */
  class EdgeS03DLocalCoefficients
#ifdef DUNE_VIRTUAL_SHAPEFUNCTIONS
    : public LocalCoefficientsInterface
#else
    : public LocalCoefficientsInterface<EdgeS03DLocalCoefficients>
#endif
  {
  public:
    //! \brief Standard constructor
    EdgeS03DLocalCoefficients () : li(6)
    {
      for (int i=0; i<6; i++)
        li[i] = LocalKey(i,1,0);
    }

    //! number of coefficients
    int size () const
    {
      return 6;
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

#endif // DUNE_EDGES03DLOCALCOEFFICIENTS_HH
