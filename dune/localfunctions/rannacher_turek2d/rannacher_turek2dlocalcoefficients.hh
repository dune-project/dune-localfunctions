// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_RANNACHER_TUREK2DLOCALCOEFFICIENTS_HH
#define DUNE_RANNACHER_TUREK2DLOCALCOEFFICIENTS_HH

#include <cstddef>
#include <vector>

#include "../common/localkey.hh"

namespace Dune {

  /** \todo Please doc me!
      \implements Dune::LocalCoefficientsVirtualImp
   */
  class RannacherTurek2DLocalCoefficients
  {
  public:
    RannacherTurek2DLocalCoefficients () : li(4)  {
      for (std::size_t i=0; i<4; i++) li[i] = LocalKey(i,1,0);
    }

    //! number of coefficients
    std::size_t size () const {
      return 4;
    }

    //! map index i to local key
    const LocalKey& localKey (std::size_t i) const {
      return li[i];
    }

  private:
    std::vector<LocalKey> li;
  };

} //namespace Dune

#endif //DUNE_RANNACHER_TUREK2DLOCALCOEFFICIENTS_HH
