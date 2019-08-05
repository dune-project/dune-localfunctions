// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LOCALFUNCTIONS_LAGRANGE_LAGRANGESIMPLEX_HH
#define DUNE_LOCALFUNCTIONS_LAGRANGE_LAGRANGESIMPLEX_HH

#include <array>
#include <numeric>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/math.hh>

#include <dune/localfunctions/common/localbasis.hh>
#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include <dune/localfunctions/common/localinterpolation.hh>
#include <dune/localfunctions/common/localkey.hh>

namespace Dune { namespace Impl
{
   /** \brief Lagrange shape functions of arbitrary order on the reference simplex

     Lagrange shape functions of arbitrary order have the property that
     \f$\hat\phi^i(x_j) = \delta_{i,j}\f$ for certain points \f$x_j\f$.

     \tparam D Type to represent the field in the domain
     \tparam R Type to represent the field in the range
     \tparam dim Dimension of the domain simplex
     \tparam k Polynomial order
   */
  template<class D, class R, unsigned int dim, unsigned int k>
  class LagrangeSimplexLocalBasis
  {
  public:
    using Traits = LocalBasisTraits<D,dim,FieldVector<D,dim>,R,1,FieldVector<R,1>,FieldMatrix<R,1,dim> >;

    /** \brief Number of shape functions
     *
     * See https://en.wikipedia.org/wiki/Figurate_number for an explanation of the formula
     */
    static constexpr unsigned int size ()
    {
      return binomial(k+dim,dim);
    }

    //! \brief Evaluate all shape functions
    void evaluateFunction(const typename Traits::DomainType& x,
                          std::vector<typename Traits::RangeType>& out) const
    {
      out.resize(size());

      // Specialization for zero-order case
      if (k==0)
      {
        out[0] = 1;
        return;
      }

      // Specialization for first-order case
      if (k==1)
      {
        out[0] = 1.0;
        for (size_t i=0; i<dim; i++)
        {
          out[0]  -= x[i];
          out[i+1] = x[i];
        }
        return;
      }

      if (dim==1)
      {
        assert(k>=2);

        auto lagrangeNode = [](unsigned int i) { return ((D)i)/k; };

        for (unsigned int i=0; i<size(); i++)
        {
          out[i] = 1.0;
          for (unsigned int alpha=0; alpha<i; alpha++)
            out[i] *= (x[0]-lagrangeNode(alpha))/(lagrangeNode(i)-lagrangeNode(alpha));
          for (unsigned int gamma=i+1; gamma<=k; gamma++)
            out[i] *= (x[0]-lagrangeNode(gamma))/(lagrangeNode(i)-lagrangeNode(gamma));
        }
        return;
      }

      if (dim!=3)
        DUNE_THROW(NotImplemented, "LagrangeSimplexLocalBasis for k>=2 only implemented for dim==1 or dim==3");

      typename Traits::DomainType kx = x;
      kx *= k;
      unsigned int n = 0;
      unsigned int i[4];
      R factor[4];
      for (i[2] = 0; i[2] <= k; ++i[2])
      {
        factor[2] = 1.0;
        for (unsigned int j = 0; j < i[2]; ++j)
          factor[2] *= (kx[2]-j) / (i[2]-j);
        for (i[1] = 0; i[1] <= k - i[2]; ++i[1])
        {
          factor[1] = 1.0;
          for (unsigned int j = 0; j < i[1]; ++j)
            factor[1] *= (kx[1]-j) / (i[1]-j);
          for (i[0] = 0; i[0] <= k - i[1] - i[2]; ++i[0])
          {
            factor[0] = 1.0;
            for (unsigned int j = 0; j < i[0]; ++j)
              factor[0] *= (kx[0]-j) / (i[0]-j);
            i[3] = k - i[0] - i[1] - i[2];
            D kx3 = k - kx[0] - kx[1] - kx[2];
            factor[3] = 1.0;
            for (unsigned int j = 0; j < i[3]; ++j)
              factor[3] *= (kx3-j) / (i[3]-j);
            out[n++] = factor[0] * factor[1] * factor[2] * factor[3];
          }
        }
      }
    }

    /** \brief Evaluate Jacobian of all shape functions
     *
     * \param x Point in the reference simplex where to evaluation the Jacobians
     * \param[out] out The Jacobians of all shape functions at the point x
     */
    void evaluateJacobian(const typename Traits::DomainType& x,
                          std::vector<typename Traits::JacobianType>& out) const
    {
      out.resize(size());

      // Specialization for k==0
      if (k==0)
      {
        std::fill(out[0][0].begin(), out[0][0].end(), 0);
        return;
      }

      // Specialization for k==1
      if (k==1)
      {
        std::fill(out[0][0].begin(), out[0][0].end(), -1);

        for (unsigned int i=0; i<dim; i++)
          for (unsigned int j=0; j<dim; j++)
            out[i+1][0][j] = (i==j);

        return;
      }

      // Specialization for dim==1
      if (dim==1)
      {
        auto lagrangeNode = [](unsigned int i) { return ((D)i)/k; };

        for (unsigned int i=0; i<=k; i++)
        {
          // x_0 derivative
          out[i][0][0] = 0.0;
          R factor=1.0;
          for (unsigned int a=0; a<i; a++)
          {
            R product=factor;
            for (unsigned int alpha=0; alpha<i; alpha++)
              product *= (alpha==a) ? 1.0/(lagrangeNode(i)-lagrangeNode(alpha))
                                    : (x[0]-lagrangeNode(alpha))/(lagrangeNode(i)-lagrangeNode(alpha));
            for (unsigned int gamma=i+1; gamma<=k; gamma++)
              product *= (lagrangeNode(gamma)-x[0])/(lagrangeNode(gamma)-lagrangeNode(i));
            out[i][0][0] += product;
          }
          for (unsigned int c=i+1; c<=k; c++)
          {
            R product=factor;
            for (unsigned int alpha=0; alpha<i; alpha++)
              product *= (x[0]-lagrangeNode(alpha))/(lagrangeNode(i)-lagrangeNode(alpha));
            for (unsigned int gamma=i+1; gamma<=k; gamma++)
              product *= (gamma==c) ? -1.0/(lagrangeNode(gamma)-lagrangeNode(i))
                                    : (lagrangeNode(gamma)-x[0])/(lagrangeNode(gamma)-lagrangeNode(i));
            out[i][0][0] += product;
          }
        }
        return;
      }

      if (dim!=3)
        DUNE_THROW(NotImplemented, "LagrangeSimplexLocalBasis only implemented for dim==3!");

      // Specialization for arbitrary order and dim==3
      typename Traits::DomainType kx = x;
      kx *= k;
      unsigned int n = 0;
      unsigned int i[4];
      R factor[4];
      for (i[2] = 0; i[2] <= k; ++i[2])
      {
        factor[2] = 1.0;
        for (unsigned int j = 0; j < i[2]; ++j)
          factor[2] *= (kx[2]-j) / (i[2]-j);
        for (i[1] = 0; i[1] <= k - i[2]; ++i[1])
        {
          factor[1] = 1.0;
          for (unsigned int j = 0; j < i[1]; ++j)
            factor[1] *= (kx[1]-j) / (i[1]-j);
          for (i[0] = 0; i[0] <= k - i[1] - i[2]; ++i[0])
          {
            factor[0] = 1.0;
            for (unsigned int j = 0; j < i[0]; ++j)
              factor[0] *= (kx[0]-j) / (i[0]-j);
            i[3] = k - i[0] - i[1] - i[2];
            D kx3 = k - kx[0] - kx[1] - kx[2];
            R sum3 = 0.0;
            factor[3] = 1.0;
            for (unsigned int j = 0; j < i[3]; ++j)
              factor[3] /= i[3] - j;
            R prod_all = factor[0] * factor[1] * factor[2] * factor[3];
            for (unsigned int j = 0; j < i[3]; ++j)
            {
              R prod = prod_all;
              for (unsigned int l = 0; l < i[3]; ++l)
                if (j == l)
                  prod *= -R(k);
                else
                  prod *= kx3 - l;
              sum3 += prod;
            }
            for (unsigned int j = 0; j < i[3]; ++j)
              factor[3] *= kx3 - j;
            for (unsigned int m = 0; m < 3; ++m)
            {
              out[n][0][m] = sum3;
              for (unsigned int j = 0; j < i[m]; ++j)
              {
                R prod = factor[3];
                for (unsigned int p = 0; p < 3; ++p)
                {
                  if (m == p)
                    for (unsigned int l = 0; l < i[p]; ++l)
                      prod *= (j==l) ? R(k) / (i[p]-l) : (kx[p]-l) / (i[p]-l);
                  else
                    prod *= factor[p];
                }
                out[n][0][m] += prod;
              }
            }
            n++;
          }
        }
      }
    }

    /** \brief Evaluate partial derivatives of any order of all shape functions
     *
     * \param order Order of the partial derivatives, in the classic multi-index notation
     * \param in Position where to evaluate the derivatives
     * \param[out] out The desired partial derivatives
     */
    void partial(const std::array<unsigned int,dim>& order,
                 const typename Traits::DomainType& in,
                 std::vector<typename Traits::RangeType>& out) const
    {
      auto totalOrder = std::accumulate(order.begin(), order.end(), 0);

      out.resize(size());

      if (totalOrder == 0) {
        evaluateFunction(in, out);
        return;
      }

      if (k==0)
      {
        out[0] = 0;
        return;
      }

      if (k==1)
      {
        if (totalOrder==1)
        {
          auto direction = std::find(order.begin(), order.end(), 1);

          out[0] = -1;
          for (unsigned int i=0; i<dim; i++)
            out[i+1] = (i==(direction-order.begin()));
        }
        else  // all higher order derivatives are zero
          std::fill(out.begin(), out.end(), 0);
        return;
      }

      DUNE_THROW(NotImplemented, "Desired derivative order is not implemented");
    }

    //! \brief Polynomial order of the shape functions
    static constexpr unsigned int order ()
    {
      return k;
    }
  };

  /** \brief Associations of the Lagrange degrees of freedom to subentities of the reference simplex
   *
   * \tparam dim Dimension of the reference simplex
   * \tparam k Polynomial order of the Lagrange space
   */
  template<unsigned int dim, unsigned int k>
  class LagrangeSimplexLocalCoefficients
  {
  public:
    //! \brief Default constructor
    LagrangeSimplexLocalCoefficients ()
    : localKeys_(size())
    {
      if (k==0)
      {
        localKeys_[0] = LocalKey(0,0,0);
        return;
      }

      if (k==1)
      {
        for (std::size_t i=0; i<size(); i++)
          localKeys_[i] = LocalKey(i,dim,0);
        return;
      }

      if (dim==1)
      {
        // Order is at least 2 here
        localKeys_[0] = LocalKey(0,1,0);          // vertex dof
        for (unsigned int i=1; i<k; i++)
          localKeys_[i] = LocalKey(0,0,i-1);      // element dofs
        localKeys_.back() = LocalKey(1,1,0);      // vertex dof
        return;
      }

      if (dim==3)
      {
        std::array<unsigned int, dim+1> vertexMap;
        for (unsigned int i=0; i<=dim; i++)
          vertexMap[i] = i;
        generateLocalKeys(vertexMap);
        return;
      }
      DUNE_THROW(NotImplemented, "LagrangeSimplexLocalCoefficients only implemented for k<=1 or dim==1 or dim==3!");
    }

    /** Constructor for variants with permuted vertices
     *
     * \param vertexmap The permutation of the vertices.  This
     *   can for instance be generated from the global indices of
     *   the vertices by reducing those to the integers 0...dim
     */
    LagrangeSimplexLocalCoefficients (const std::array<unsigned int, dim+1> vertexMap)
    : localKeys_(size())
    {
      if (dim!=3)
        DUNE_THROW(NotImplemented, "LagrangeSimplexLocalCoefficients only implemented for dim==3!");

      generateLocalKeys(vertexMap);
    }

    //! number of coefficients
    static constexpr std::size_t size ()
    {
      return binomial(k+dim,dim);
    }

    //! get i'th index
    const LocalKey& localKey (std::size_t i) const
    {
      return localKeys_[i];
    }

  private:
    std::vector<LocalKey> localKeys_;

    void generateLocalKeys(const std::array<unsigned int, dim+1> vertexMap)
    {
      if (k==0)
      {
        localKeys_[0] = LocalKey(0,0,0);
        return;
      }
      if (dim!=3)
        DUNE_THROW(NotImplemented, "LagrangeSimplexLocalCoefficients only implemented for dim==3!");

      unsigned int subindex[16];
      unsigned int codim_count[4] = {0};
      for (unsigned int m = 1; m < 16; ++m)
      {
        unsigned int codim = !(m&1) + !(m&2) + !(m&4) + !(m&8);
        subindex[m] = codim_count[codim]++;
      }

      int a1 = (3*k + 12)*k + 11;
      int a2 = -3*k - 6;
      unsigned int dof_count[16] = {0};
      unsigned int i[4];
      for (i[3] = 0; i[3] <= k; ++i[3])
        for (i[2] = 0; i[2] <= k - i[3]; ++i[2])
          for (i[1] = 0; i[1] <= k - i[2] - i[3]; ++i[1])
          {
            i[0] = k - i[1] - i[2] - i[3];
            unsigned int j[4];
            unsigned int entity = 0;
            unsigned int codim = 0;
            for (unsigned int m = 0; m < 4; ++m)
            {
              j[m] = i[vertexMap[m]];
              entity += !!j[m] << m;
              codim += !j[m];
            }
            int local_index = j[3]*(a1 + (a2 + j[3])*j[3])/6
                              + j[2]*(2*(k - j[3]) + 3 - j[2])/2 + j[1];
            localKeys_[local_index] = LocalKey(subindex[entity], codim, dof_count[entity]++);
          }
    }
  };

  /** \brief Evaluate the degrees of freedom of a Lagrange basis
   *
   * \tparam LocalBasis The corresponding set of shape functions
   */
  template<class LocalBasis>
  class LagrangeSimplexLocalInterpolation
  {
    static const int kdiv = (LocalBasis::order() == 0 ? 1 : LocalBasis::order());
  public:

    /** \brief Evaluate a given function at the Lagrange nodes
     *
     * \tparam F Type of function to evaluate
     * \tparam C Type used for the values of the function
     * \param[in] ff Function to evaluate
     * \param[out] out Array of function values
     */
    template<typename F, typename C>
    void interpolate (const F& ff, std::vector<C>& out) const
    {
      constexpr auto dim = LocalBasis::Traits::dimDomain;
      constexpr auto k = LocalBasis::order();
      using D = typename LocalBasis::Traits::DomainFieldType;

      typename LocalBasis::Traits::DomainType x;
      auto&& f = Impl::makeFunctionWithCallOperator<typename LocalBasis::Traits::DomainType>(ff);

      out.resize(LocalBasis::size());

      // Specialization for zero-order case
      if (k==0)
      {
        auto center = ReferenceElements<D,dim>::simplex().position(0,0);
        out[0] = f(center);
        return;
      }

      // Specialization for first-order case
      if (k==1)
      {
        // vertex 0
        std::fill(x.begin(), x.end(), 0);
        out[0] = f(x);

        // remaining vertices
        for (int i=0; i<dim; i++)
        {
          for (int j=0; j<dim; j++)
            x[j] = (i==j);

          out[i+1] = f(x);
        }
        return;
      }

      if (dim==1)
      {
        for (unsigned int i=0; i<k+1; i++)
        {
          x[0] = ((D)i)/k;
          out[i] = f(x);
        }
        return;
      }

      if (dim!=3)
        DUNE_THROW(NotImplemented, "LagrangeSimplexLocalInterpolation only implemented for dim==1 or dim==3!");

      int n=0;
      for (int i2 = 0; i2 <= (int)k; i2++)
        for (int i1 = 0; i1 <= (int)k-i2; i1++)
          for (int i0 = 0; i0 <= (int)k-i1-i2; i0++)
          {
            x[0] = ((D)i0)/((D)kdiv);
            x[1] = ((D)i1)/((D)kdiv);
            x[2] = ((D)i2)/((D)kdiv);
            out[n] = f(x);
            n++;
          }
    }

  };

} }    // namespace Dune::Impl

namespace Dune
{
  /** \brief Lagrange finite element for simplices with arbitrary compile-time dimension and polynomial order
   *
   * \tparam D type used for domain coordinates
   * \tparam R type used for function values
   * \tparam d dimension of the reference element
   * \tparam k polynomial order
   */
  template<class D, class R, int d, int k>
  class LagrangeSimplexLocalFiniteElement
  {
  public:
    /** \brief Export number types, dimensions, etc.
     */
    using Traits = LocalFiniteElementTraits<Impl::LagrangeSimplexLocalBasis<D,R,d,k>,
                                            Impl::LagrangeSimplexLocalCoefficients<d,k>,
                                            Impl::LagrangeSimplexLocalInterpolation<Impl::LagrangeSimplexLocalBasis<D,R,d,k> > >;

    /** \brief Returns the local basis, i.e., the set of shape functions
     */
    const typename Traits::LocalBasisType& localBasis () const
    {
      return basis_;
    }

    /** \brief Returns the assignment of the degrees of freedom to the element subentities
     */
    const typename Traits::LocalCoefficientsType& localCoefficients () const
    {
      return coefficients_;
    }

    /** \brief Returns object that evaluates degrees of freedom
     */
    const typename Traits::LocalInterpolationType& localInterpolation () const
    {
      return interpolation_;
    }

    /** \brief The number of shape functions */
    static constexpr std::size_t size ()
    {
      return Impl::LagrangeSimplexLocalBasis<D,R,d,k>::size();
    }

    /** \brief The reference element that the local finite element is defined on
     */
    static constexpr GeometryType type ()
    {
      return GeometryTypes::simplex(d);
    }

  private:
    Impl::LagrangeSimplexLocalBasis<D,R,d,k> basis_;
    Impl::LagrangeSimplexLocalCoefficients<d,k> coefficients_;
    Impl::LagrangeSimplexLocalInterpolation<Impl::LagrangeSimplexLocalBasis<D,R,d,k> > interpolation_;
  };

}        // namespace Dune

#endif   // DUNE_LOCALFUNCTIONS_LAGRANGE_LAGRANGESIMPLEX_HH
