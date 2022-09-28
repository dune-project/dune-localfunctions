// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_BREZZIDOUGLASFORTINMARINI_CUBE_LOCALINTERPOLATION_HH
#define DUNE_LOCALFUNCTIONS_BREZZIDOUGLASFORTINMARINI_CUBE_LOCALINTERPOLATION_HH

#include <algorithm>
#include <array>
#include <bitset>
#include <vector>
#include <limits>

#include <dune/common/fvector.hh>
#include <dune/common/math.hh>
#include <dune/common/rangeutilities.hh>
#include <dune/common/typetraits.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/localfunctions/common/localinterpolation.hh>


namespace Dune
{

  /**
   * \ingroup LocalInterpolationImplementation
   * \brief Interpolation for Brezzi-Douglas-Fortin-Marini shape functions on cubes.
   *
   * \tparam D      Type of represent the field in the domain.
   * \tparam R      Type of represent the field in the domain.
   * \tparam dim    dimension of the reference element, must be >= 2.
   * \tparam order  order of the element, must be >= 1.
   *
   * \nosubgrouping
   */
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

    /**
     * \brief compute the i'th shifted Legendre function on [0,1]
     *
     * \param i  index of function to compute
     * \param x  Position
     */
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

    /**
     * \brief calculate the i`th multi index in graded lexigraphic order
     *
     * \tparam d     number of components of the multi index
     * \tparam kMax  maximum absolute value to consider. (default is `unlimited`)
     *
     * \param i  index
     *
     */
    template<std::size_t d, std::size_t kMax = std::numeric_limits<std::size_t>::max()>
    constexpr inline static auto unrank (std::size_t i)
      -> std::array<std::size_t, d>
    {
      assert( i < binomial(d+kMax, kMax));

      std::array<std::size_t, d> mi{};
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

    /**
     * \brief embed the a local position on a face into the reference element
     *
     * \param face  index of the face
     * \param x     local position
     */
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
    //! \brief Standard constructor
    BDFMCubeLocalInterpolation ()
    {
      std::fill(sign_.begin(), sign_.end(), 1.0);
    }

    /**
     * \brief Make set number s, where 0 <= s < 2^(2*dim)
     *
     * \param s  Edge orientation indicator
     */
    BDFMCubeLocalInterpolation (std::bitset<numFaces> s)
    {
      for (auto i : range(numFaces))
      {
        sign_[i] = s[i] ? -1 : 1;
        normal_[i][i/2] = i%2 ? 1 : -1;
      }
    }

    /**
     * \brief Interpolate a given function with shape functions
     *
     * \tparam F  Function type for function which should be interpolated
     * \tparam C  Coefficient vector type
     *
     * \param ff   function which should be interpolated
     * \param out  return value, vector of coefficients
     */
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

    /**
     * \brief Interpolate a given function with shape functions on a face of the reference element
     *
     * \tparam F  Function type for function which should be interpolated
     * \tparam C  Coefficient vector type
     *
     * \param face  index of the face on which to interpolate
     * \param f     function which should be interpolated
     * \param out   return value, vector of coefficients
     *
     */
    template<class F, class C>
    void trace (unsigned int face, const F& f, C& out) const
    {
      assert(out.size() >= numFaces*faceDofs + interiorDofs);
      assert(face < numFaces);

      const auto o = face*faceDofs;
      const auto& n = normal_[face];
      const auto& s = sign_[face];
      const auto& c = n[face/2];

      const auto& rule = QuadratureRules<DomainFieldType, dim-1>::rule(GeometryTypes::cube(dim-1), order+(order-1)+5);
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

          out[o+i] += (y*n) * phi * qp.weight() * (i%2 ? c : s);
        }
      }
    }

    /**
     * \brief Interpolate a given function with shape functions in the interior of the reference element
     *
     * \tparam F  Function type for function which should be interpolated
     * \tparam C  Coefficient vector type
     *
     * \param ff   function which should be interpolated
     * \param out  return value, vector of coefficients
     */
    template<class F, class C>
    void interior (const F& f, C& out) const
    {
      assert(out.size() >= numFaces*faceDofs + interiorDofs);

       const auto o = numFaces*faceDofs;

      const auto& rule = QuadratureRules<DomainFieldType, dim>::rule(GeometryTypes::cube(dim), order+std::max((int)order-2,(int)0)+5);
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


  template<class D, class R, unsigned int dim, unsigned int order>
  constexpr std::size_t BDFMCubeLocalInterpolation<D, R, dim, order>::numFaces;

  template<class D, class R, unsigned int dim, unsigned int order>
  constexpr unsigned int BDFMCubeLocalInterpolation<D, R, dim, order>::interiorDofs;

  template<class D, class R, unsigned int dim, unsigned int order>
  constexpr unsigned int BDFMCubeLocalInterpolation<D, R, dim, order>::faceDofs;


#ifndef DOXYGEN
  template<class D, class R, unsigned int dim>
  class BDFMCubeLocalInterpolation<D, R, dim, 0>
  {
    static_assert(AlwaysFalse<D>::value,
                  "`BDFMCubeLocalCoefficients` not defined for order 0.");
  };
#endif //#ifndef DOXYGEN

} // namespace Dune

#endif // #ifndef DUNE_LOCALFUNCTIONS_BREZZIDOUGLASFORTINMARINI_CUBE_LOCALINTERPOLATION_HH
