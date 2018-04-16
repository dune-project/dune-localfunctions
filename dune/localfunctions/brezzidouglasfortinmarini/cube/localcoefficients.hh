#ifndef DUNE_LOCALFUNCTIONS_BREZZIDOUGLASFORTINMARINI_CUBE_LOCALCOEFFICIENTS_HH
#define DUNE_LOCALFUNCTIONS_BREZZIDOUGLASFORTINMARINI_CUBE_LOCALCOEFFICIENTS_HH

#include <cstddef>
#include <vector>

#include <dune/common/mathutilities.hh>
#include <dune/common/rangeutilities.hh>
#include <dune/common/typetraits.hh>

#include <dune/localfunctions/common/localkey.hh>

namespace Dune
{

  template<class D, class R, unsigned int dim, unsigned int order>
  class BDFMCubeLocalCoefficients
  {
    static constexpr unsigned int interiorDofs = dim*binomial(dim+order-2, order-2);
    static constexpr unsigned int faceDofs     = binomial(dim+order-2, order-1);

    static constexpr std::size_t numFaces = 2*dim;
    static constexpr std::size_t numDofs  = numFaces*faceDofs + interiorDofs;

    public:
      BDFMCubeLocalCoefficients () : li(numDofs)
      {
        for (auto j : range(numFaces))
          for (auto i : range(faceDofs))
            li[j*faceDofs + i] = LocalKey(j, 1, i);

        for (auto i : range(interiorDofs))
          li[numFaces*faceDofs + i] = LocalKey(0, 0, i);
      }

      std::size_t size () const { return numDofs; }
      auto localKey (std::size_t i) const -> const LocalKey& { return li[i]; }

    private:
      std::vector<LocalKey> li;
  };


  template<class D, class R, unsigned int dim>
  class BDFMCubeLocalCoefficients<D, R, dim, 0>
  {
    static_assert( AlwaysFalse<D>::value,
                   "`BDFMCubeLocalCoefficients` not defined for order 0." );
  };

} // namespace Dune

#endif // #ifndef DUNE_LOCALFUNCTIONS_BREZZIDOUGLASFORTINMARINI_CUBE_LOCALCOEFFICIENTS_HH
