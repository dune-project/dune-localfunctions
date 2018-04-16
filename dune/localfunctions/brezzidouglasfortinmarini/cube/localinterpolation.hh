#ifndef DUNE_LOCALFUNCTIONS_BREZZIDOUGLASFORTINMARINI_CUBE_LOCALINTERPOLATION_HH
#define DUNE_LOCALFUNCTIONS_BREZZIDOUGLASFORTINMARINI_CUBE_LOCALINTERPOLATION_HH

#include <algorithm>
#include <array>
#include <bitset>
#include <vector>

#include <dune/common/fvector.hh>
#include <dune/common/mathutilities.hh>
#include <dune/common/rangeutilities.hh>
#include <dune/common/typetraits.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/localfunctions/common/localinterpolation.hh>


namespace Dune
{

  template<class D, class R, unsigned int dim, unsigned int order>
  class BDFMCubeLocalInterpolation
  {
    using DomainType      = FieldVector<D, dim>;
    using FaceDomainType  = FieldVector<D, dim-1>;
    using RangeType       = FieldVector<R, dim>;
    using DomainFieldType = D;
    using RangeFieldType  = R;

    static constexpr std::size_t numFaces = 2*dim;

    static constexpr unsigned int interiorDofs = dim*binomial(dim+order-2, order-2);
    static constexpr unsigned int faceDofs     = binomial(dim+order-2, order-1);

    inline static auto legendre( unsigned int i, const DomainFieldType& x )
      -> RangeFieldType
    {
      switch( i )
      {
        case 0:
          return 1;
        case 1:
          return 2*x-1;
        case 2:
          return 6*x*x-6*x+1;
        case 3:
          return 20*x*x*x-30*x*x+12*x-1;
        default:
          return ((2*i-1)*(2*x-1)*legendre(x,i-1) - (i-1)*legendre(x,i-2)) / i;
      }
    }

    template<std::size_t k>
    using MultiIndex = FieldVector<std::size_t, k>;

    template<std::size_t d, std::size_t kMax = -1>
    constexpr inline static auto unrank (std::size_t i)
      -> MultiIndex<d>
    {
      assert( i < binomial(d+kMax, kMax));

      MultiIndex<d> mi;
      if (i == 0)
        return mi;

      std::size_t k = 0, b = 1;
      for(;k <= kMax && b <= i; ++k, b = binomial(d+k-1, k))
        i -= b;

      std::size_t p = 0;
      for(; p < d && i > 0; ++p)
      {
        std::size_t m = 0, c = 1;
        for(; m <= k && c <= i; ++m, c = binomial(d-p+m-2, m))
          i -= c;

        mi[p] = k-m;
        k = m;
      }
      mi[p] = k;

      return mi;
    }

    inline static auto embed (unsigned int face, const FaceDomainType& x )
      -> DomainType
    {
      DomainType y;

      std::size_t j = 0;
      for (auto i : range(dim))
        if (i == face/2)
          y[i] = face%2;
        else
          y[i] = x[j++];

      return y;
    }

  public:
    BDFMCubeLocalInterpolation ()
    {
      std::fill(sign_.begin(), sign_.end(), 1.0);
    }

    BDFMCubeLocalInterpolation (std::bitset<numFaces> s)
    {
      for (auto i : range(numFaces))
      {
        sign_[i] = s[i] ? -1 : 1;
        normal_[i][i/2] = i%2 ? 1 : -1;
      }
    }

    template<class F, class C>
    void interpolate (const F& ff, C& out) const
    {
      out.resize(numFaces*faceDofs + interiorDofs);
      std::fill(out.begin(),out.end(), 0.0);

      auto&& f = Impl::makeFunctionWithCallOperator<DomainType>(ff);

      for(auto i : range(numFaces))
        trace(i, f, out);

      interior(f, out);
    }

    template<class F, class C>
    void trace (unsigned int face, const F& f, C& out) const
    {
      assert(out.size() >= numFaces*faceDofs + interiorDofs);
      assert(face < numFaces);

      const auto o = face*faceDofs;
      const auto& n = normal_[face];
      const auto& s = sign_[face];

      const auto& rule = QuadratureRules<DomainFieldType, dim-1>::rule(GeometryTypes::cube(dim-1), order+1);
      for (const auto& qp : rule)
      {
        const auto& x = qp.position();
        auto y = f(embed(face, x));

        for (auto i : range(faceDofs))
        {
          auto mi = unrank<dim-1,order-1>(i);

          RangeFieldType phi = 1.;
          for (auto j : range(dim-1))
            phi *= legendre(mi[j], x[j]);

          out[o+i] += (y*n) * phi * qp.weight() * (i%2 ? R(1) : s);
        }
      }
    }

    template<class F, class C>
    void interior (const F& f, C& out) const
    {
      assert(out.size() >= numFaces*faceDofs + interiorDofs);

       const auto o = numFaces*faceDofs;

      const auto& rule = QuadratureRules<DomainFieldType, dim>::rule(GeometryTypes::cube(dim), order+1);
      for(const auto& qp : rule)
      {
        const auto& x = qp.position();
        auto y = f(x);

        for (auto i : range(interiorDofs/dim))
        {
          auto mi = unrank<dim,order-2>(i);

          RangeFieldType phi = 1.;
          for (auto j : range(dim))
            phi *= legendre(mi[j], x[j]);

          for (auto j : range(dim))
            out[o+dim*i+j] +=  y[j]* phi * qp.weight();
        }
      }
    }

  private:
    std::array<RangeFieldType, numFaces> sign_;
    std::array<DomainType, numFaces> normal_;
  };

  template<class D, class R, unsigned int dim>
  class BDFMCubeLocalInterpolation<D, R, dim, 0>
  {
    static_assert(AlwaysFalse<D>::value,
                  "`BDFMCubeLocalCoefficients` not defined for order 0.");
  };

} // namespace Dune

#endif // #ifndef DUNE_LOCALFUNCTIONS_BREZZIDOUGLASFORTINMARINI_CUBE_LOCALINTERPOLATION_HH
