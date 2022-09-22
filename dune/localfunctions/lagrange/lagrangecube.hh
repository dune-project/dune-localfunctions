// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
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
  // Forward declaration
  template<class LocalBasis>
  class LagrangeCubeLocalInterpolation;

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
    friend class LagrangeCubeLocalInterpolation<LagrangeCubeLocalBasis<D,R,dim,k> >;

    // i-th Lagrange polynomial of degree k in one dimension
    static R p(unsigned int i, D x)
    {
      R result(1.0);
      for (unsigned int j=0; j<=k; j++)
        if (j!=i) result *= (k*x-j)/((int)i-(int)j);
      return result;
    }

    // derivative of ith Lagrange polynomial of degree k in one dimension
    static R dp(unsigned int i, D x)
    {
      R result(0.0);

      for (unsigned int j=0; j<=k; j++)
      {
        if (j!=i)
        {
          R prod( (k*1.0)/((int)i-(int)j) );
          for (unsigned int l=0; l<=k; l++)
            if (l!=i && l!=j)
              prod *= (k*x-l)/((int)i-(int)l);
          result += prod;
        }
      }
      return result;
    }

    // Second derivative of j-th Lagrange polynomial of degree k in one dimension
    // Formula and notation taken from https://en.wikipedia.org/wiki/Lagrange_polynomial#Derivatives
    static R ddp(unsigned int j, D x)
    {
      R result(0.0);

      for (unsigned int i=0; i<=k; i++)
      {
        if (i==j)
          continue;

        R sum(0);

        for (unsigned int m=0; m<=k; m++)
        {
          if (m==i || m==j)
            continue;

          R prod( (k*1.0)/((int)j-(int)m) );
          for (unsigned int l=0; l<=k; l++)
            if (l!=i && l!=j && l!=m)
              prod *= (k*x-l)/((int)j-(int)l);
          sum += prod;
        }

        result += sum * ( (k*1.0)/((int)j-(int)i) );
      }

      return result;
    }

    // Return i as a d-digit number in the (k+1)-nary system
    static std::array<unsigned int,dim> multiindex (unsigned int i)
    {
      std::array<unsigned int,dim> alpha;
      for (unsigned int j=0; j<dim; j++)
      {
        alpha[j] = i % (k+1);
        i = i/(k+1);
      }
      return alpha;
    }

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

      // General case
      for (size_t i=0; i<size(); i++)
      {
        // convert index i to multiindex
        std::array<unsigned int,dim> alpha(multiindex(i));

        // initialize product
        out[i] = 1.0;

        // dimension by dimension
        for (unsigned int j=0; j<dim; j++)
          out[i] *= p(alpha[j],x[j]);
      }
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

      // The general case

      // Loop over all shape functions
      for (size_t i=0; i<size(); i++)
      {
        // convert index i to multiindex
        std::array<unsigned int,dim> alpha(multiindex(i));

        // Loop over all coordinate directions
        for (unsigned int j=0; j<dim; j++)
        {
          // Initialize: the overall expression is a product
          // if j-th bit of i is set to -1, else 1
          out[i][0][j] = dp(alpha[j],x[j]);

          // rest of the product
          for (unsigned int l=0; l<dim; l++)
            if (l!=j)
              out[i][0][j] *= p(alpha[l],x[l]);
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
        else if (totalOrder == 2)
        {

          for (size_t i=0; i<size(); i++)
          {
            // convert index i to multiindex
            std::array<unsigned int,dim> alpha(multiindex(i));

            // Initialize: the overall expression is a product
            out[i][0] = 1.0;

            // rest of the product
            for (std::size_t l=0; l<dim; l++)
            {
              switch (order[l])
              {
                case 0:
                  out[i][0] *= p(alpha[l],in[l]);
                  break;
                case 1:
                  //std::cout << "dp: " << dp(alpha[l],in[l]) << std::endl;
                  out[i][0] *= dp(alpha[l],in[l]);
                  break;
                case 2:
                  out[i][0] *= 0;
                  break;
                default:
                  DUNE_THROW(NotImplemented, "Desired derivative order is not implemented");
              }
            }
          }
        }
        else
          DUNE_THROW(NotImplemented, "Partial derivative of order " << totalOrder << " is not implemented!");

        return;
      }

      // The case k>1

      // Loop over all shape functions
      for (size_t i=0; i<size(); i++)
      {
        // convert index i to multiindex
        std::array<unsigned int,dim> alpha(multiindex(i));

        // Initialize: the overall expression is a product
        out[i][0] = 1.0;

        // rest of the product
        for (std::size_t l=0; l<dim; l++)
        {
          switch (order[l])
          {
            case 0:
              out[i][0] *= p(alpha[l],in[l]);
              break;
            case 1:
              out[i][0] *= dp(alpha[l],in[l]);
              break;
            case 2:
              out[i][0] *= ddp(alpha[l],in[l]);
              break;
            default:
              DUNE_THROW(NotImplemented, "Desired derivative order is not implemented");
          }
        }
      }
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
    // Return i as a d-digit number in the (k+1)-nary system
    static std::array<unsigned int,dim> multiindex (unsigned int i)
    {
      std::array<unsigned int,dim> alpha;
      for (unsigned int j=0; j<dim; j++)
      {
        alpha[j] = i % (k+1);
        i = i/(k+1);
      }
      return alpha;
    }

    /** \brief Set the 'subentity' field for each dof for a 1d element */
    void setup1d(std::vector<unsigned int>& subEntity)
    {
      assert(k>0);

      unsigned lastIndex=0;

      /* edge and vertex numbering
         0----0----1
       */

      // edge (0)
      subEntity[lastIndex++] = 0;                 // corner 0
      for (unsigned i = 0; i < k - 1; ++i)
        subEntity[lastIndex++] = 0;               // inner dofs of element (0)

      subEntity[lastIndex++] = 1;                 // corner 1

      assert(power(k+1,dim)==lastIndex);
    }

    void setup2d(std::vector<unsigned int>& subEntity)
    {
      assert(k>0);

      unsigned lastIndex=0;

      // LocalKey: entity number, entity codim, dof indices within each entity
      /* edge and vertex numbering
       2----3----3
       |         |
       |         |
       0         1
       |         |
       |         |
       0----2----1
       */

      // lower edge (2)
      subEntity[lastIndex++] = 0;                 // corner 0
      for (unsigned i = 0; i < k - 1; ++i)
        subEntity[lastIndex++] = 2;           // inner dofs of lower edge (2)

      subEntity[lastIndex++] = 1;                 // corner 1

      // iterate from bottom to top over inner edge dofs
      for (unsigned e = 0; e < k - 1; ++e) {
        subEntity[lastIndex++] = 0;                   // left edge (0)
        for (unsigned i = 0; i < k - 1; ++i)
          subEntity[lastIndex++] = 0;                     // face dofs
        subEntity[lastIndex++] = 1;                   // right edge (1)
      }

      // upper edge (3)
      subEntity[lastIndex++] = 2;                 // corner 2
      for (unsigned i = 0; i < k - 1; ++i)
        subEntity[lastIndex++] = 3;                   // inner dofs of upper edge (3)

      subEntity[lastIndex++] = 3;                 // corner 3

      assert(power(k+1,dim)==lastIndex);
    }

    void setup3d(std::vector<unsigned int>& subEntity)
    {
      assert(k>0);

      unsigned lastIndex=0;
#ifndef NDEBUG
      const unsigned numIndices = power(k+1,dim);
      const unsigned numFaceIndices = power(k+1,dim-1);
#endif
      const unsigned numInnerEdgeDofs = k-1;

      // LocalKey: entity number, entity codim, dof indices within each entity
      /* edge and vertex numbering

              6---(11)--7              6---------7
             /|        /|             /|  (5)   /|
           (8)|      (9)|            / | top   / |
           / (2)     / (3)          /  |(3)bac/k |
          4---(10)--5   |          4---------5   |
          |   |     |   |      left|(0)|     |(1)|right
          |   2--(7)|---3          |   2-----|---3
         (0) /     (1) /           |(2)front |  /
          |(4)      |(5)           | /  (4)  | /
          |/        |/             |/ bottom |/
          0---(6)---1              0---------1
       */

      // bottom face (4)
      lastIndex=0;
      // lower edge (6)
      subEntity[lastIndex++] = 0;              // corner 0
      for (unsigned i = 0; i < numInnerEdgeDofs; ++i)
        subEntity[lastIndex++] = 6;                // inner dofs of lower edge (6)

      subEntity[lastIndex++] = 1;              // corner 1

      // iterate from bottom to top over inner edge dofs
      for (unsigned e = 0; e < numInnerEdgeDofs; ++e) {
        subEntity[lastIndex++] = 4;                // left edge (4)
        for (unsigned i = 0; i < numInnerEdgeDofs; ++i)
          subEntity[lastIndex++] = 4;                       // inner face dofs
        subEntity[lastIndex++] = 5;                 // right edge (5)
      }

      // upper edge (7)
      subEntity[lastIndex++] = 2;              // corner 2
      for (unsigned i = 0; i < k - 1; ++i)
        subEntity[lastIndex++] = 7;                // inner dofs of upper edge (7)
      subEntity[lastIndex++] = 3;                // corner 3

      assert(numFaceIndices==lastIndex);       // added 1 face so far
      /////////////////////////////////////////// end bottom face (4)

      ///////////////////// inner faces
      for(unsigned f = 0; f < numInnerEdgeDofs; ++f) {

        // lower edge (connecting  edges 0 and 1)
        subEntity[lastIndex++] = 0;                // dof on edge 0
        for (unsigned i = 0; i < numInnerEdgeDofs; ++i)
          subEntity[lastIndex++] = 2;                            // dof in front face
        subEntity[lastIndex++] = 1;                // dof on edge 1

        // iterate from bottom to top over inner edge dofs
        for (unsigned e = 0; e < numInnerEdgeDofs; ++e) {
          subEntity[lastIndex++] = 0;                  // on left face (0)
          for (unsigned i = 0; i < numInnerEdgeDofs; ++i)
            subEntity[lastIndex++] = 0;                    // volume dofs
          subEntity[lastIndex++] = 1;                  // right face (1)
        }

        // upper edge (connecting  edges 0 and 1)
        subEntity[lastIndex++] = 2;                // dof on edge 2
        for (unsigned i = 0; i < numInnerEdgeDofs; ++i)
          subEntity[lastIndex++] = 3;                  // dof on rear face (3)
        subEntity[lastIndex++] = 3;                // dof on edge 3

        assert(lastIndex==(f+1+1)*numFaceIndices);
      }

      ////////////////////////////////////////// top face (5)
      // lower edge (10)
      subEntity[lastIndex++] = 4;              // corner 4
      for (unsigned i = 0; i < k - 1; ++i)
        subEntity[lastIndex++] = 10;                // inner dofs on lower edge (10)
      subEntity[lastIndex++] = 5;              // corner 5

      // iterate from bottom to top over inner edge dofs
      for (unsigned e = 0; e < k - 1; ++e) {
        subEntity[lastIndex++] = 8;                // left edge (8)
        for (unsigned i = 0; i < k - 1; ++i)
          subEntity[lastIndex++] = 5;                  // face dofs
        subEntity[lastIndex++] = 9;                // right edge (9)
      }

      // upper edge (11)
      subEntity[lastIndex++] = 6;              // corner 6
      for (unsigned i = 0; i < k - 1; ++i)
        subEntity[lastIndex++] = 11;                // inner dofs of upper edge (11)
      subEntity[lastIndex++] = 7;              // corner 7

      assert(numIndices==lastIndex);
    }

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

      // Now: the general case

      // Set up array of codimension-per-dof-number
      std::vector<unsigned int> codim(size());

      for (std::size_t i=0; i<codim.size(); i++)
      {
        codim[i] = 0;

        // Codimension gets increased by 1 for each coordinate direction
        // where dof is on boundary
        std::array<unsigned int,dim> mIdx = multiindex(i);
        for (unsigned int j=0; j<dim; j++)
          if (mIdx[j]==0 or mIdx[j]==k)
            codim[i]++;
      }

      // Set up index vector (the index of the dof in the set of dofs of a given subentity)
      // Algorithm: the 'index' has the same ordering as the dof number 'i'.
      // To make it consecutive we interpret 'i' in the (k+1)-adic system, omit all digits
      // that correspond to axes where the dof is on the element boundary, and transform the
      // rest to the (k-1)-adic system.
      std::vector<unsigned int> index(size());

      for (std::size_t i=0; i<size(); i++)
      {
        index[i] = 0;

        std::array<unsigned int,dim> mIdx = multiindex(i);

        for (int j=dim-1; j>=0; j--)
          if (mIdx[j]>0 && mIdx[j]<k)
            index[i] = (k-1)*index[i] + (mIdx[j]-1);
      }

      // Set up entity and dof numbers for each (supported) dimension separately
      std::vector<unsigned int> subEntity(size());

      if (dim==1) {

        setup1d(subEntity);

      } else if (dim==2) {

        setup2d(subEntity);

      } else if (dim==3) {

        setup3d(subEntity);

      } else
        DUNE_THROW(Dune::NotImplemented, "LagrangeCubeLocalCoefficients for order " << k << " and dim == " << dim);

      for (size_t i=0; i<size(); i++)
        localKeys_[i] = LocalKey(subEntity[i], codim[i], index[i]);
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

      // The general case
      for (unsigned int i=0; i<LocalBasis::size(); i++)
      {
        // convert index i to multiindex
        std::array<unsigned int,dim> alpha(LocalBasis::multiindex(i));

        // Generate coordinate of the i-th Lagrange point
        for (unsigned int j=0; j<dim; j++)
          x[j] = (1.0*alpha[j])/k;

        out[i] = f(x);
      }
    }

  };

} }    // namespace Dune::Impl

namespace Dune
{
  /** \brief Lagrange finite element for cubes with arbitrary compile-time dimension and polynomial order
   *
   * \tparam D Type used for domain coordinates
   * \tparam R Type used for function values
   * \tparam dim dimension of the reference element
   * \tparam k Polynomial order in one coordinate direction
   */
  template<class D, class R, int dim, int k>
  class LagrangeCubeLocalFiniteElement
  {
  public:
    /** \brief Export number types, dimensions, etc.
     */
    using Traits = LocalFiniteElementTraits<Impl::LagrangeCubeLocalBasis<D,R,dim,k>,
                                            Impl::LagrangeCubeLocalCoefficients<dim,k>,
                                            Impl::LagrangeCubeLocalInterpolation<Impl::LagrangeCubeLocalBasis<D,R,dim,k> > >;

    /** \brief Default constructor
     *
     * \deprecated This explicit implementation only exists to work around a bug in clang 3.8
     *   which disappeared in clang 6
     */
    LagrangeCubeLocalFiniteElement() {}

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
      return power(k+1,dim);
    }

    /** \brief The reference element that the local finite element is defined on
     */
    static constexpr GeometryType type ()
    {
      return GeometryTypes::cube(dim);
    }

  private:
    Impl::LagrangeCubeLocalBasis<D,R,dim,k> basis_;
    Impl::LagrangeCubeLocalCoefficients<dim,k> coefficients_;
    Impl::LagrangeCubeLocalInterpolation<Impl::LagrangeCubeLocalBasis<D,R,dim,k> > interpolation_;
  };

}        // namespace Dune

#endif   // DUNE_LOCALFUNCTIONS_LAGRANGE_LAGRANGECUBE_HH
