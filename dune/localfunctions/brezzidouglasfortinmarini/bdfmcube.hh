// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LOCALFUNCTIONS_BREZZIDOUGLASFORTINMARINI_BDFMCUBE_HH
#define DUNE_LOCALFUNCTIONS_BREZZIDOUGLASFORTINMARINI_BDFMCUBE_HH

#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/localfiniteelementtraits.hh>

#include <dune/localfunctions/brezzidouglasfortinmarini/cube/localbasis.hh>
#include <dune/localfunctions/brezzidouglasfortinmarini/cube/localcoefficients.hh>
#include <dune/localfunctions/brezzidouglasfortinmarini/cube/localinterpolation.hh>


namespace Dune
{

  /**
   * \brief Brezzi-Douglas-Fortin-Marini finite elements for cubes
   *
   * \ingroup BrezziDouglasFortinMarini
   *
   * \tparam D      Type to represent the field in the domain.
   * \tparam R      Type to represent the field in the range.
   * \tparam dim    dimension of the reference elements, must be >= 2.
   * \tparam order  order of the element, must be >= 1.
   */
  template<class D, class R, unsigned int dim, unsigned int order>
  class BDFMCubeLocalFiniteElement
  {
    using LocalBasis          = BDFMCubeLocalBasis<D, R, dim, order>;
    using LocalCoefficients   = BDFMCubeLocalCoefficients<D, R, dim, order>;
    using LocalInterpolation  = BDFMCubeLocalInterpolation<D, R, dim, order>;

  public:
    using Traits = LocalFiniteElementTraits<LocalBasis, LocalCoefficients, LocalInterpolation  >;

    //! \brief Standard constructor
    BDFMCubeLocalFiniteElement () {}

    /**
     * \brief Make set number  s, where 0 <= s < 2^(2*dim)
     *
     * \param s  Edge orientation indicator
     */
    BDFMCubeLocalFiniteElement (std::bitset<2*dim> s)
      : basis( s ), interpolation( s )
    {}

    auto localBasis () const -> const LocalBasis& { return basis; }
    auto localCoefficients () const -> const LocalCoefficients& { return coefficients; }
    auto localInterpolation () const -> const LocalInterpolation& { return interpolation; }

    /** \brief Number of shape functions in this finite element */
    unsigned int size () const { return basis.size(); }
    static constexpr auto type () -> GeometryType { return GeometryTypes::cube(dim); }

  private:
    LocalBasis basis;
    LocalCoefficients coefficients;
    LocalInterpolation interpolation;
  };

} // namespace Dune

#endif // #ifndef DUNE_LOCALFUNCTIONS_BREZZIDOUGLASFORTINMARINI_BDFMCUBE_HH
