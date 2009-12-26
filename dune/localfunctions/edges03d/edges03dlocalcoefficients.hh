// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_EDGES03DLOCALCOEFFICIENTS_HH
#define DUNE_EDGES03DLOCALCOEFFICIENTS_HH

#include <cstddef>
#include <iostream>
#include <vector>

#include "../common/localcoefficients.hh"

namespace Dune
{

  /**@ingroup LocalLayoutImplementation
     \brief Layout map for lowest order edge elements on tetrahedrons

     \nosubgrouping
     \implements Dune::LocalCoefficientsVirtualImp
   */
  class EdgeS03DLocalCoefficients
#ifdef DUNE_VIRTUAL_SHAPEFUNCTIONS
    : public LocalCoefficientsInterface
#endif
  {
  public:
    //! \brief Standard constructor
    EdgeS03DLocalCoefficients () : li(6)
    {
      for (int i=0; i<6; i++)
        li[i] = LocalKey(i,2,0);
    }

    //! number of coefficients
    std::size_t size () const
    {
      return 6;
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

#endif // DUNE_EDGES03DLOCALCOEFFICIENTS_HH
