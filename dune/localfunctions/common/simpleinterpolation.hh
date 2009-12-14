// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_SIMPLEINTERPOLATION_HH
#define DUNE_SIMPLEINTERPOLATION_HH

#include <vector>

#include "interpolation.hh"

namespace Dune
{

  //! Wrap a LocalInterpolation class and turn it into a "global" Interpolation
  /**
   *  \tparam LI Type of the LocalInterpolation
   */
  template<class LI>
  class SimpleInterpolation
    : public InterpolationInterface<SimpleInterpolation<LI> >
  {
    //! store a reference to the wrapped LocalInterpolation object.
    const LI& li;

  public:
    //! Construct a SimpleInterpolation object
    /**
     *  \param [in] li_ A reference to the LocalInterpolation object.  This
     *                  reference is stored internally and should be valid for
     *                  as long as the interpolate() method of this
     *                  SimpleInterpolation object is used.
     */
    SimpleInterpolation(const LI& li_)
      : li(li_)
    {}

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
      li.interpolate(f, out);
    }
  };

}

#endif // DUNE_SIMPLEINTERPOLATION_HH
