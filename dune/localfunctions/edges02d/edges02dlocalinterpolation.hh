// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_EDGES02DLOCALINTERPOLATION_HH
#define DUNE_EDGES02DLOCALINTERPOLATION_HH

#include <cmath>

#include <dune/grid/common/genericreferenceelements.hh>

namespace Dune
{
  /**@ingroup LocalBasisImplementation
         \brief Interpolation for experimental lowest order edge elements for triangles.

         \tparam LB LocalBasisImplementation

     \note This class does not implement the usual LocalInterpolationInterface
           since that does not make much sense for vector valued elements.  An
           experimental interface providing global rather than local values is
           provided instead.  Be warned that this interface is subject to
           change without notice, however.

         \nosubgrouping
   */

  template<class LB>
  class EdgeS02DLocalInterpolation
  {
  public:
    //! contruct an interpolation instance with default orientations
    EdgeS02DLocalInterpolation()
    {
      s[0] = 1; s[1] = 1; s[2] = 1;
    }

    //! contruct an interpolation instance with the given orientations
    //! \param orientations Bit-map of orientations for each shape function;
    //! bit 0 = 0 means default orientation for the first shape function, bit
    //! 0 = 1 means inverted orientation for the first shape function.
    EdgeS02DLocalInterpolation(unsigned int orientations)
    {
      s[0] = 1; s[1] = 1; s[2] = 1;
      for(int i = 0; i < 3; ++i)
        if(orientations & (1<<i)) s[i] *= -1;
    }

    //! \brief Local interpolation of a function
    /**
     * \tparam F        Type of function to interpolate.  The class should
     *                  provide a method <tt>void evaluate(const DomainType
     *                  &x, RangeType &y) const</tt> which is used to evaluate
     *                  the function on the reference element.  This method
     *                  should expect local coordinates on the reference
     *                  element as input <tt>x</tt> but should return global
     *                  values in <tt>y</tt>.
     * \tparam C        Type of coefficients.
     * \tparam Geometry Type of Geometry.
     *
     * \param[in]  f        Function instance used to interpolate.
     * \param[out] out      Resulting coefficients vector.
     * \param[in]  geometry geometry of the element used for interpolation.
     */
    template<typename F, typename C, typename Geometry>
    void interpolateGlobal (const F& f, std::vector<C>& out,
                            const Geometry &geometry) const
    {
      static const GenericReferenceElement<typename LB::Traits::DomainFieldType, 2> &refElem
        = GenericReferenceElements<typename LB::Traits::DomainFieldType, 2>::simplex();

      typename LB::Traits::DomainType vertex[3];
      for(int i = 0; i < 3; ++i)
        vertex[i] = geometry.corner(i);

      typename LB::Traits::DomainType tangent;
      typename LB::Traits::DomainType x;
      typename LB::Traits::RangeType y;

      out.resize(3);
      for(int j = 0; j < 3; ++j) {
        int v0 = refElem.subEntity(j, 1, 0, 2);
        int v1 = refElem.subEntity(j, 1, 1, 2);
        if(v0 > v1) std::swap(v0, v1);

        tangent = vertex[v1]; tangent -= vertex[v0];
        x = vertex[v1]; x += vertex[v0]; x /= 2;

        f.evaluate(x,y);

        out[j] = s[j] * (tangent * y) / tangent.two_norm();
      }
    }

  private:
    //! The signs
    typename LB::Traits::RangeFieldType s[3];
  };
}

#endif // DUNE_EDGES02DLOCALINTERPOLATION_HH
