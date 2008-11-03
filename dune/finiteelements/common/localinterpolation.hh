// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LOCALINTERPOLATION_HH
#define DUNE_LOCALINTERPOLATION_HH

#include <iostream>
#include <vector>

namespace Dune
{

  template<class Imp>
  class LocalInterpolationInterface
  {
  public:

    //! determine coefficients interpolating a given function
    template<typename F, typename C>
    void interpolate (const F& f, std::vector<C>& out) const
    {
      asImp().interpolate(f,out);
    }

  private:
    Imp& asImp () {return static_cast<Imp &> (*this);}
    const Imp& asImp () const {return static_cast<const Imp &>(*this);}
  };

}

#endif
