// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LOCALINTERPOLATION_HH
#define DUNE_LOCALINTERPOLATION_HH

#include <vector>

namespace Dune
{

  //! Interface class for interpolating a local basis
  template<class Imp>
  class LocalInterpolationInterface
  {
  public:

    //! determine coefficients interpolating a given function
    /**
     * \tparam F Type of function to interpolate.  The class should provide a
     *           method <tt>void evaluate(const DomainType &x, RangeType &y)
     *           const</tt> which is used to evaluate the function on the
     *           reference element.
     * \tparam C Type of coefficients.
     *
     * \param[in]  f   Function instance used to interpolate.
     * \param[out] out Resulting coefficients vector.
     */
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
