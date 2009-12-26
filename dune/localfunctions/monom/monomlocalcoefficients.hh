// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_MONOMLOCALCOEFFICIENTS_HH
#define DUNE_MONOMLOCALCOEFFICIENTS_HH

#include <cstddef>
#include <iostream>
#include <vector>

#include "../common/localcoefficients.hh"

namespace Dune
{

  /**@ingroup LocalLayoutImplementation
         \brief Layout map for Monom elements

         \nosubgrouping
   */
  template<int static_size>
  class MonomLocalCoefficients
#ifdef DUNE_VIRTUAL_SHAPEFUNCTIONS
    : public LocalCoefficientsInterface
#endif
  {
  public:
    //! \brief Standard constructor
    MonomLocalCoefficients ()
      : index(static_size, LocalKey(0,0,0))
    {
      for(int i = 0; i < static_size; ++i)
        index[i].index(i);
    }

    //! number of coefficients
    std::size_t size () const
    {
      return static_size;
    }

    //! get i'th index
    const LocalKey& localKey (std::size_t i) const
    {
      return index[i];
    }

  private:
    std::vector<LocalKey> index;
  };

}
#endif //DUNE_MONOMLOCALCOEFFICIENTS_HH
