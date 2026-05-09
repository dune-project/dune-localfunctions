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
#include <dune/common/typetraits.hh>

#include <dune/geometry/referenceelements.hh>

#include <dune/localfunctions/common/localbasis.hh>
#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include <dune/localfunctions/common/localkey.hh>

namespace Dune { namespace Impl
{
  // The traits provide static or dynamic order and size information
  template<unsigned int dim, int compileTimeOrder>
  struct LagrangeCubeOrderTraits;

  // The traits provide static order and size information
  template<unsigned int dim, int compileTimeOrder>
  requires (compileTimeOrder >= 0)
  struct LagrangeCubeOrderTraits<dim,compileTimeOrder>
  {
    static constexpr bool is_static_order = true;

    constexpr LagrangeCubeOrderTraits(int /*runTimeOrder*/ = compileTimeOrder) {}

    /**
     * \brief Number of shape functions
     */
    static constexpr std::size_t size ()
    {
      return power(compileTimeOrder+1, dim);
    }

    /**
     * \brief Polynomial order of the one-dimensional shape functions
     */
    static constexpr unsigned int order ()
    {
      return compileTimeOrder;
    }

    // Return i as a d-digit number in the (k+1)-nary system
    static constexpr std::array<unsigned int,dim> multiindex (unsigned int i)
    {
      std::array<unsigned int,dim> alpha;
      for (unsigned int j=0; j<dim; j++)
      {
        alpha[j] = i % (order()+1);
        i = i/(order()+1);
      }
      return alpha;
    }
  };


  template<unsigned int dim>
  struct LagrangeCubeOrderTraits<dim,-1>
  {
    std::size_t size_;
    unsigned int order_;

    static constexpr bool is_static_order = false;

    /**
     * \brief Constructor computes the size and stores the order
     */
    constexpr LagrangeCubeOrderTraits (int runTimeOrder)
      : size_(power(runTimeOrder+1, dim))
      , order_(runTimeOrder)
    {
      if (runTimeOrder < 0)
        DUNE_THROW(Dune::InvalidStateException, "LagrangeCube: run-time order must be >= 0");
    }

    /**
     * \brief Number of shape functions
     */
    constexpr std::size_t size () const
    {
      return size_;
    }

    /**
     * \brief Polynomial order of the one-dimensional shape functions
     */
    constexpr unsigned int order () const
    {
      return order_;
    }

    // Return i as a d-digit number in the (k+1)-nary system
    constexpr std::array<unsigned int,dim> multiindex (unsigned int i) const
    {
      std::array<unsigned int,dim> alpha;
      for (unsigned int j=0; j<dim; j++)
      {
        alpha[j] = i % (order()+1);
        i = i/(order()+1);
      }
      return alpha;
    }
  };

   /** \brief Lagrange shape functions of arbitrary order on the reference cube [0,1]^d

     Lagrange shape functions of arbitrary order have the property that
     \f$\hat\phi^i(x_j) = \delta_{i,j}\f$ for certain points \f$x_j\f$.

     \tparam D Type to represent the field in the domain
     \tparam R Type to represent the field in the range
     \tparam dim Dimension of the domain cube
     \tparam compileTimeOrder Polynomial order; -1 means "order provided at run-time"
   */
  template<class D, class R, unsigned int dim, int compileTimeOrder>
  class LagrangeCubeLocalBasis
    : private LagrangeCubeOrderTraits<dim,compileTimeOrder>
  {
    using OrderTraits = LagrangeCubeOrderTraits<dim,compileTimeOrder>;

    // i-th Lagrange polynomial of degree k in one dimension
    constexpr R p (unsigned int i, D x) const
    {
      R result(1.0);
      for (unsigned int j=0; j<=order(); j++)
        if (j!=i) result *= (order()*x-j)/((int)i-(int)j);
      return result;
    }

    // derivative of ith Lagrange polynomial of degree k in one dimension
    constexpr R dp (unsigned int i, D x) const
    {
      R result(0.0);

      for (unsigned int j=0; j<=order(); j++)
      {
        if (j!=i)
        {
          R prod( (order()*1.0)/((int)i-(int)j) );
          for (unsigned int l=0; l<=order(); l++)
            if (l!=i && l!=j)
              prod *= (order()*x-l)/((int)i-(int)l);
          result += prod;
        }
      }
      return result;
    }

    // Second derivative of j-th Lagrange polynomial of degree k in one dimension
    // Formula and notation taken from https://en.wikipedia.org/wiki/Lagrange_polynomial#Derivatives
    constexpr R ddp(unsigned int j, D x) const
    {
      R result(0.0);

      for (unsigned int i=0; i<=order(); i++)
      {
        if (i==j)
          continue;

        R sum(0);

        for (unsigned int m=0; m<=order(); m++)
        {
          if (m==i || m==j)
            continue;

          R prod( (order()*1.0)/((int)j-(int)m) );
          for (unsigned int l=0; l<=order(); l++)
            if (l!=i && l!=j && l!=m)
              prod *= (order()*x-l)/((int)j-(int)l);
          sum += prod;
        }

        result += sum * ( (order()*1.0)/((int)j-(int)i) );
      }

      return result;
    }

    using OrderTraits::multiindex;

  public:
    using Traits = LocalBasisTraits<D,dim,FieldVector<D,dim>,R,1,FieldVector<R,1>,FieldMatrix<R,1,dim> >;

    constexpr LagrangeCubeLocalBasis () requires (OrderTraits::is_static_order)
      : OrderTraits()
    {}

    /** \brief Constructor from OrderTraits
     */
    explicit constexpr LagrangeCubeLocalBasis (OrderTraits orderTraits)
      : OrderTraits(orderTraits)
    {}

    using OrderTraits::size;
    using OrderTraits::order;

    //! \brief Evaluate all shape functions
    constexpr void evaluateFunction (const typename Traits::DomainType& x,
                                     std::vector<typename Traits::RangeType>& out) const
    {
      out.resize(size());

      // Specialization for zero-order case
      if (order()==0)
      {
        out[0] = 1;
        return;
      }

      if (order()==1)
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
    constexpr void evaluateJacobian (const typename Traits::DomainType& x,
                                     std::vector<typename Traits::JacobianType>& out) const
    {
      out.resize(size());

      // Specialization for common case
      if (order()==0)
      {
        std::fill(out[0][0].begin(), out[0][0].end(), 0);
        return;
      }

      // Specialization for common case
      if (order()==1)
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
     * \param partialOrders Orders of the partial derivatives, in the classic multi-index notation
     * \param in Position where to evaluate the derivatives
     * \param[out] out The desired partial derivatives
     */
    constexpr void partial (const std::array<unsigned int,dim>& partialOrders,
                            const typename Traits::DomainType& in,
                            std::vector<typename Traits::RangeType>& out) const
    {
      auto totalOrder = std::accumulate(partialOrders.begin(), partialOrders.end(), 0);

      out.resize(size());

      if (order()==0)
      {
        out[0] = (totalOrder==0);
        return;
      }

      if (order()==1)
      {
        if (totalOrder == 0)
        {
          evaluateFunction(in, out);
        }
        else if (totalOrder == 1)
        {
          out.resize(size());

          auto direction = std::distance(partialOrders.begin(),
                                         std::find(partialOrders.begin(), partialOrders.end(), 1));
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
              switch (partialOrders[l])
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
          switch (partialOrders[l])
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
  };

  /** \brief Associations of the Lagrange degrees of freedom to subentities of the reference cube
   *
   * \tparam dim Dimension of the reference cube
   * \tparam compileTimeOrder Polynomial order of the Lagrange space in one direction
   */
  template<unsigned int dim, int compileTimeOrder>
  class LagrangeCubeLocalCoefficients
    : private LagrangeCubeOrderTraits<dim,compileTimeOrder>
  {
    using OrderTraits = LagrangeCubeOrderTraits<dim,compileTimeOrder>;

    using OrderTraits::multiindex;

    /** \brief Set the 'subentity' field for each dof for a 1d element */
    void setup1d (std::vector<unsigned int>& subEntity)
    {
      const unsigned int k = order();
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

    void setup2d (std::vector<unsigned int>& subEntity)
    {
      const unsigned int k = order();
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

    void setup3d (std::vector<unsigned int>& subEntity)
    {
      const unsigned int k = order();
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

    /** \brief Initializes the localKeys_ array
     */
    void setup ()
    {
      const unsigned int k = order();
      localKeys_.resize(size());

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

  public:
    //! \brief Default constructor
    explicit LagrangeCubeLocalCoefficients (OrderTraits orderTraits)
      : OrderTraits(orderTraits)
    {
      setup();
    }

    using OrderTraits::size;
    using OrderTraits::order;

    //! get i-th index
    constexpr const LocalKey& localKey (std::size_t i) const
    {
      return localKeys_[i];
    }

  private:
    std::vector<LocalKey> localKeys_;
  };

  /** \brief Evaluate the degrees of freedom of a Lagrange basis
   *
   * \tparam D Type to represent the field in the domain
   * \tparam R Type to represent the field in the range
   * \tparam dim Dimension of the domain cube
   * \tparam compileTimeOrder Polynomial order; -1 means "order provided at run-time"
   */
  template<class D, class R, unsigned int dim, int compileTimeOrder>
  class LagrangeCubeLocalInterpolation
    : private LagrangeCubeOrderTraits<dim, compileTimeOrder>
  {
    using OrderTraits = LagrangeCubeOrderTraits<dim, compileTimeOrder>;
    using Traits = typename LagrangeCubeLocalBasis<D,R,dim,compileTimeOrder>::Traits;

    using OrderTraits::multiindex;
    using OrderTraits::size;
    using OrderTraits::order;

  public:
    /** \brief Constructor for a given set of shape functions */
    explicit constexpr LagrangeCubeLocalInterpolation (OrderTraits orderTraits)
      : OrderTraits(orderTraits)
    {}

    /** \brief Evaluate a given function at the Lagrange nodes
     *
     * \tparam F Type of function to evaluate
     * \tparam C Type used for the values of the function
     * \param[in] f Function to evaluate
     * \param[out] out Array of function values
     */
    template<typename F, typename C>
    constexpr void interpolate (const F& f, std::vector<C>& out) const
    {
      const unsigned int k = order();
      out.resize(size());

      // Specialization for zero-order case
      if (k==0)
      {
        auto center = ReferenceElements<D,dim>::cube().position(0,0);
        out[0] = f(center);
        return;
      }

      typename Traits::DomainType x;

      // Specialization for first-order case
      if (k==1)
      {
        for (unsigned int i=0; i<size(); i++)
        {
          // Generate coordinate of the i-th corner of the reference cube
          for (unsigned int j=0; j<dim; j++)
            x[j] = (i & (1<<j)) ? 1.0 : 0.0;

          out[i] = f(x);
        }
        return;
      }

      // The general case
      for (unsigned int i=0; i<size(); i++)
      {
        // convert index i to multiindex
        std::array<unsigned int,dim> alpha(multiindex(i));

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
   * \tparam compileTimeOrder Polynomial order in one coordinate direction.
   *           The default -1 means "order provided at run-time"
   */
  template<class D, class R, int dim, int compileTimeOrder = -1>
  class LagrangeCubeLocalFiniteElement
    : private Impl::LagrangeCubeOrderTraits<dim, compileTimeOrder>
  {
    using OrderTraits = Impl::LagrangeCubeOrderTraits<dim, compileTimeOrder>;

  public:
    /** \brief Export number types, dimensions, etc.
     */
    using Traits = LocalFiniteElementTraits<
      Impl::LagrangeCubeLocalBasis<D,R,dim,compileTimeOrder>,
      Impl::LagrangeCubeLocalCoefficients<dim,compileTimeOrder>,
      Impl::LagrangeCubeLocalInterpolation<D,R,dim,compileTimeOrder> >;

    //! \brief Constructor for static order
    constexpr LagrangeCubeLocalFiniteElement ()
      : OrderTraits()
      , basis_(*this)
      , coefficients_(*this)
      , interpolation_(*this)
    {
      static_assert(OrderTraits::is_static_order, "Default constructor only allowed for compile-time order >= 0");
      static_assert(compileTimeOrder >= 0, "Default constructor only allowed for compile-time order >= 0");
    }

    //! \brief Constructor for dynamic order
    explicit constexpr LagrangeCubeLocalFiniteElement (int runTimeOrder)
      : OrderTraits(runTimeOrder)
      , basis_(*this)
      , coefficients_(*this)
      , interpolation_(*this)
    {
      if (runTimeOrder < 0)
        DUNE_THROW(Dune::InvalidStateException, "LagrangeCubeLocalFiniteElement: run-time order must be non-negative!");
      if (compileTimeOrder >= 0 && compileTimeOrder != runTimeOrder)
        DUNE_THROW(Dune::InvalidStateException, "LagrangeCubeLocalFiniteElement: Compile-time order must be identical to run-time order!");
    }

    /** \brief Returns the local basis, i.e., the set of shape functions
     */
    constexpr const typename Traits::LocalBasisType& localBasis () const
    {
      return basis_;
    }

    /** \brief Returns the assignment of the degrees of freedom to the element subentities
     */
    constexpr const typename Traits::LocalCoefficientsType& localCoefficients () const
    {
      return coefficients_;
    }

    /** \brief Returns object that evaluates degrees of freedom
     */
    constexpr const typename Traits::LocalInterpolationType& localInterpolation () const
    {
      return interpolation_;
    }

    /** \brief The number of shape functions */
    using OrderTraits::size;

    /** \brief The reference element that the local finite element is defined on
     */
    static constexpr GeometryType type ()
    {
      return GeometryTypes::cube(dim);
    }

  private:
    typename Traits::LocalBasisType basis_;
    typename Traits::LocalCoefficientsType coefficients_;
    typename Traits::LocalInterpolationType interpolation_;
  };

}        // namespace Dune

#endif   // DUNE_LOCALFUNCTIONS_LAGRANGE_LAGRANGECUBE_HH
