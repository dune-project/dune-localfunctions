// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_REFINED_P0_LOCALCOEFFICIENTS_HH
#define DUNE_REFINED_P0_LOCALCOEFFICIENTS_HH

#include <cstddef>
#include <iostream>
#include <vector>

#include "../common/localcoefficients.hh"

namespace Dune
{
  template<int base, int exponent>
  struct Pow
  {
    enum {value = base * Pow<base,exponent-1>::value};
  };

  template<int base>
  struct Pow<base,0>
  {
    enum {value = 1};
  };

  /**@ingroup LocalLayoutImplementation
     \brief Layout map for RefinedP0 elements

     \nosubgrouping
     \implements Dune::LocalCoefficientsVirtualImp
   */
  template<unsigned int k>
  class RefinedP0LocalCoefficients
#if DUNE_VIRTUAL_SHAPEFUNCTIONS
    : public LocalCoefficientsInterface
#endif
  {
    enum {N = Pow<2,k>::value};

  public:
    RefinedP0LocalCoefficients () :
      localKeys_(N)
    {
      // All functions are associated to the element
      for (int i = 0; i < N; ++i)
        localKeys_[i] = LocalKey(0,0,i);
    }

    //! number of coefficients
    std::size_t size () const
    {
      return N;
    }

    //! get i'th index
    const LocalKey& localKey (std::size_t i) const
    {
      return localKeys_[i];
    }

  private:
    std::vector<LocalKey> localKeys_;

  };

}

#endif
