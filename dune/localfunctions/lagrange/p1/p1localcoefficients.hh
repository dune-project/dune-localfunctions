// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_P1_LOCALCOEFFICIENTS_HH
#define DUNE_P1_LOCALCOEFFICIENTS_HH

#include <cstddef>
#include <iostream>
#include <vector>

#include <dune/localfunctions/common/localkey.hh>

namespace Dune
{

  /** \brief Local coefficients for simplex P1 elements

     \nosubgrouping
     \implements Dune::LocalCoefficientsVirtualImp
   */
  template <int dim>
  class P1LocalCoefficients
  {
  public:
    //! \brief Standard constructor
    P1LocalCoefficients () : li(size())
    {
      for (std::size_t i=0; i<size(); i++)
        li[i] = LocalKey(i,dim,0);
    }

    //! number of coefficients
    std::size_t size () const
    {
      return dim+1;
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
