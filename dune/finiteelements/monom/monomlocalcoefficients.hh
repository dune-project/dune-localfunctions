// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_MONOMLOCALCOEFFICIENTS_HH
#define DUNE_MONOMLOCALCOEFFICIENTS_HH

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
  class MonomLocalCoefficients : public LocalCoefficientsInterface<MonomLocalCoefficients<static_size> >
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
    int size () const
    {
      return static_size;
    }

    //! get i'th index
    const LocalKey& localKey (int i) const
    {
      return index[i];
    }

  private:
    std::vector<LocalKey> index;
  };

}
#endif //DUNE_MONOMLOCALCOEFFICIENTS_HH
