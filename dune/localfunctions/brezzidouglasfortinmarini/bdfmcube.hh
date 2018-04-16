#ifndef DUNE_LOCALFUNCTIONS_BREZZIDOUGLASFORTINMARINI_BDFMCUBE_HH
#define DUNE_LOCALFUNCTIONS_BREZZIDOUGLASFORTINMARINI_BDFMCUBE_HH

#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/localfiniteelementtraits.hh>

#include <dune/localfunctions/brezzidouglasfortinmarini/cube/localbasis.hh>
#include <dune/localfunctions/brezzidouglasfortinmarini/cube/localcoefficients.hh>
#include <dune/localfunctions/brezzidouglasfortinmarini/cube/localinterpolation.hh>


namespace Dune
{

  template<class D, class R, unsigned int dim, unsigned int order>
  class BDFMCubeLocalFiniteElement
  {
    using LocalBasis          = BDFMCubeLocalBasis<D, R, dim, order>;
    using LocalCoefficients   = BDFMCubeLocalCoefficients<D, R, dim, order>;
    using LocalInterpolation  = BDFMCubeLocalInterpolation<D, R, dim, order>;

  public:
    using Traits = LocalFiniteElementTraits<LocalBasis, LocalCoefficients, LocalInterpolation  >;

    BDFMCubeLocalFiniteElement () {}

    BDFMCubeLocalFiniteElement (std::bitset<2*dim> s)
      : basis( s ), interpolation( s )
    {}

    auto localBasis () const -> const LocalBasis& { return basis; }
    auto localCoefficients () const -> const LocalCoefficients& { return coefficients; }
    auto localInterpolation () const -> const LocalInterpolation& { return interpolation; }

    unsigned int size () const { return basis.size(); }
    static constexpr auto type () -> GeometryType { return GeometryTypes::cube(dim); }

  private:
    LocalBasis basis;
    LocalCoefficients coefficients;
    LocalInterpolation interpolation;
  };

} // namespace Dune

#endif // #ifndef DUNE_LOCALFUNCTIONS_BREZZIDOUGLASFORTINMARINI_BDFMCUBE_HH
