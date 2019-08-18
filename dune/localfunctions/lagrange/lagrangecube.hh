// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LOCALFUNCTIONS_LAGRANGE_LAGRANGECUBE_HH
#define DUNE_LOCALFUNCTIONS_LAGRANGE_LAGRANGECUBE_HH

#include <array>
#include <numeric>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/math.hh>

#include <dune/geometry/referenceelements.hh>

#include <dune/localfunctions/common/localbasis.hh>
#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include <dune/localfunctions/common/localinterpolation.hh>
#include <dune/localfunctions/common/localkey.hh>

namespace Dune { namespace Impl
{
   /** \brief Lagrange shape functions of arbitrary order on the reference cube [0,1]^d

     Lagrange shape functions of arbitrary order have the property that
     \f$\hat\phi^i(x_j) = \delta_{i,j}\f$ for certain points \f$x_j\f$.

     \tparam D Type to represent the field in the domain
     \tparam R Type to represent the field in the range
     \tparam dim Dimension of the domain cube
     \tparam k Polynomial order
   */
  template<class D, class R, unsigned int dim, unsigned int k>
  class LagrangeCubeLocalBasis
  {
  public:
    using Traits = LocalBasisTraits<D,dim,FieldVector<D,dim>,R,1,FieldVector<R,1>,FieldMatrix<R,1,dim> >;

    /** \brief Number of shape functions
     */
    static constexpr unsigned int size ()
    {
      return power(k+1, dim);
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

      if (k==1)
      {
        for (size_t i=0; i<size(); i++)
        {
          out[i] = 1;

          for (unsigned int j=0; j<dim; j++)
            // if j-th bit of i is set multiply with x[j], else with 1-x[j]
            out[i] *= (i & (1<<j)) ? x[j] :  1-x[j];
        }
        return;
      }

      DUNE_THROW(NotImplemented, "LagrangeCubeLocalBasis");
    }

    /** \brief Evaluate Jacobian of all shape functions
     *
     * \param x Point in the reference cube where to evaluation the Jacobians
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
        // Loop over all shape functions
        for (size_t i=0; i<size(); i++)
        {
          // Loop over all coordinate directions
          for (unsigned int j=0; j<dim; j++)
          {
            // Initialize: the overall expression is a product
            // if j-th bit of i is set to 1, else -11
            out[i][0][j] = (i & (1<<j)) ? 1 : -1;

            for (unsigned int l=0; l<dim; l++)
            {
              if (j!=l)
                // if l-th bit of i is set multiply with x[l], else with 1-x[l]
                out[i][0][j] *= (i & (1<<l)) ? x[l] :  1-x[l];
            }
          }
        }
        return;
      }

      DUNE_THROW(NotImplemented, "LagrangeCubeLocalBasis");
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

      if (k==0)
      {
        out[0] = (totalOrder==0);
        return;
      }

      if (k==1)
      {
        if (totalOrder == 0)
        {
          evaluateFunction(in, out);
        }
        else if (totalOrder == 1)
        {
          out.resize(size());

          auto direction = std::distance(order.begin(), std::find(order.begin(), order.end(), 1));
          if (direction >= dim)
            DUNE_THROW(RangeError, "Direction of partial derivative not found!");

          // Loop over all shape functions
          for (std::size_t i = 0; i < size(); ++i)
          {
            // Initialize: the overall expression is a product
            // if j-th bit of i is set to 1, otherwise to -1
            out[i] = (i & (1<<direction)) ? 1 : -1;

            for (unsigned int j = 0; j < dim; ++j)
            {
              if (direction != j)
                // if j-th bit of i is set multiply with in[j], else with 1-in[j]
                out[i] *= (i & (1<<j)) ? in[j] :  1-in[j];
            }
          }
        }
        else
          DUNE_THROW(NotImplemented, "Partial derivative of order " << totalOrder << " is not implemented!");

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

  /** \brief Associations of the Lagrange degrees of freedom to subentities of the reference cube
   *
   * \tparam dim Dimension of the reference cube
   * \tparam k Polynomial order of the Lagrange space in one direction
   */
  template<unsigned int dim, unsigned int k>
  class LagrangeCubeLocalCoefficients
  {
  public:
    //! \brief Default constructor
    LagrangeCubeLocalCoefficients ()
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

      DUNE_THROW(NotImplemented, "LagrangeCubeLocalCoefficients");
    }

    //! number of coefficients
    static constexpr std::size_t size ()
    {
      return power(k+1,dim);
    }

    //! get i-th index
    const LocalKey& localKey (std::size_t i) const
    {
      return localKeys_[i];
    }

  private:
    std::vector<LocalKey> localKeys_;
  };

  /** \brief Evaluate the degrees of freedom of a Lagrange basis
   *
   * \tparam LocalBasis The corresponding set of shape functions
   */
  template<class LocalBasis>
  class LagrangeCubeLocalInterpolation
  {
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
        auto center = ReferenceElements<D,dim>::cube().position(0,0);
        out[0] = f(center);
        return;
      }

      // Specialization for first-order case
      if (k==1)
      {
        for (unsigned int i=0; i<LocalBasis::size(); i++)
        {
          // Generate coordinate of the i-th corner of the reference cube
          for (int j=0; j<dim; j++)
            x[j] = (i & (1<<j)) ? 1.0 : 0.0;

          out[i] = f(x);
        }
        return;
      }

      DUNE_THROW(NotImplemented, "LagrangeCubeLocalInterpolation");
    }

  };

} }    // namespace Dune::Impl

namespace Dune
{
  /** \brief Lagrange finite element for cubes with arbitrary compile-time dimension and polynomial order
   *
   * \tparam D Type used for domain coordinates
   * \tparam R Type used for function values
   * \tparam d dimension of the reference element
   * \tparam k Polynomial order in one coordinate direction
   */
  template<class D, class R, int d, int k>
  class LagrangeCubeLocalFiniteElement
  {
  public:
    /** \brief Export number types, dimensions, etc.
     */
    using Traits = LocalFiniteElementTraits<Impl::LagrangeCubeLocalBasis<D,R,d,k>,
                                            Impl::LagrangeCubeLocalCoefficients<d,k>,
                                            Impl::LagrangeCubeLocalInterpolation<Impl::LagrangeCubeLocalBasis<D,R,d,k> > >;

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
      return power(k+1,d);
    }

    /** \brief The reference element that the local finite element is defined on
     */
    static constexpr GeometryType type ()
    {
      return GeometryTypes::cube(d);
    }

  private:
    Impl::LagrangeCubeLocalBasis<D,R,d,k> basis_;
    Impl::LagrangeCubeLocalCoefficients<d,k> coefficients_;
    Impl::LagrangeCubeLocalInterpolation<Impl::LagrangeCubeLocalBasis<D,R,d,k> > interpolation_;
  };

}        // namespace Dune

#endif   // DUNE_LOCALFUNCTIONS_LAGRANGE_LAGRANGECUBE_HH
