// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_NEDELEC_NEDELEC1STKINDSIMPLEX_HH
#define DUNE_LOCALFUNCTIONS_NEDELEC_NEDELEC1STKINDSIMPLEX_HH

#include <numeric>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/localbasis.hh>
#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include <dune/localfunctions/common/localinterpolation.hh>   // For deprecated makeFunctionWithCallOperator
#include <dune/localfunctions/common/localkey.hh>

namespace Dune
{
namespace Impl
{
  /** \brief Nedelec shape functions of the first kind on reference simplices
   *
   * \note These shape functions are implemented for the reference simplices only!
   *   The transformation to other simplices has to be done by the user.
   *
   * \tparam D Number type used for domain coordinates
   * \tparam R Number type used for function values
   * \tparam dim Dimension of the reference element
   * \tparam k Order of the Nedelec element
   */
  template<class D, class R, int dim, int k>
  class Nedelec1stKindSimplexLocalBasis
  {
    // Number of edges of the reference simplex
    constexpr static std::size_t numberOfEdges = dim*(dim+1)/2;

  public:
    using Traits = LocalBasisTraits<D,dim,FieldVector<D,dim>,
                                    R,dim,FieldVector<R,dim>,
                                    FieldMatrix<R,dim,dim> >;

    /** \brief Constructor with default edge orientation
     *
     * Default orientation means that all edges point from the vertex with lower index
     * to the vertex with higher index. The tangential parts of all shape functions
     * point in the direction of the edge.
     */
    Nedelec1stKindSimplexLocalBasis()
    {
      std::fill(edgeOrientation_.begin(), edgeOrientation_.end(), 1.0);
    }

    /** \brief Constructor for a given edge orientation
     */
    Nedelec1stKindSimplexLocalBasis(std::bitset<numberOfEdges> edgeOrientation)
    : Nedelec1stKindSimplexLocalBasis()
    {
      for (std::size_t i=0; i<edgeOrientation_.size(); i++)
        edgeOrientation_[i] *= edgeOrientation[i] ? -1.0 : 1.0;
    }

    //! \brief Number of shape functions
    static constexpr unsigned int size()
    {
      static_assert(dim==2 || dim==3, "Nedelec shape functions are implemented only for 2d and 3d simplices.");
      if (dim==2)
        return k * (k+2);
      if (dim==3)
        return k * (k+2) * (k+3) / 2;
    }

    /** \brief Evaluate values of all shape functions at a given point
     *
     * \param[in]  in  The evaluation point
     * \param[out] out Jacobians of all shape functions at that point
     */
    void evaluateFunction(const typename Traits::DomainType& in,
                           std::vector<typename Traits::RangeType>& out) const
    {
      static_assert(k==1, "Evaluating Nédélec shape functions is implemented only for first order.");
      out.resize(size());

      if (dim==2)
      {
        // First-order Nédélec shape functions on a triangle are of the form
        //
        //         (a1, a2) + b(-x2, x1)^T,     a_1, a_2, b \in R
        out[0] = {D(1) - in[1],  in[0]};
        out[1] = {in[1],        -in[0]+D(1)};
        out[2] = {-in[1],        in[0]};
      }

      if constexpr (dim==3)
      {
        // First-order Nédélec shape functions on a tetrahedron are of the form
        //
        //          a + b \times x,       a, b \in R^3
        //
        // The following coefficients create the six basis vectors
        // that are dual to the edge degrees of freedom:
        //
        // a[0] = { 1,  0,  0}              b[0] = { 0, -1,  1}
        // a[1] = { 0,  1,  0}              b[1] = { 1,  0, -1}
        // a[2] = { 0,  0,  0}              b[2] = { 0,  0,  1}
        // a[3] = { 0,  0,  1}              b[3] = {-1,  1,  0}
        // a[4] = { 0,  0,  0}              b[4] = { 0, -1,  0}
        // a[5] = { 0,  0,  0}              b[5] = { 1,  0,  0}
        //
        // The following implementation uses these values, and simply
        // skips all the zeros.

        out[0] = { 1 - in[1] - in[2],     in[0]        ,     in[0]        };
        out[1] = {     in[1]        , 1 - in[0] - in[2],             in[1]};
        out[2] = {   - in[1]        ,     in[0]        , 0                };
        out[3] = {             in[2],             in[2], 1 - in[0] - in[1]};
        out[4] = {            -in[2], 0                ,     in[0]        };
        out[5] = { 0                ,            -in[2],             in[1]};
      }

      for (std::size_t i=0; i<out.size(); i++)
        out[i] *= edgeOrientation_[i];
    }

    /** \brief Evaluate Jacobians of all shape functions at a given point
     *
     * \param[in]  in  The evaluation point
     * \param[out] out Jacobians of all shape functions at that point
     */
    void evaluateJacobian(const typename Traits::DomainType& in,
                          std::vector<typename Traits::JacobianType>& out) const
    {
      out.resize(size());
      if (dim==2)
      {
        out[0][0] = { 0, -1};
        out[0][1] = { 1,  0};

        out[1][0] = { 0,  1};
        out[1][1] = {-1,  0};

        out[2][0] = { 0, -1};
        out[2][1] = { 1,  0};
      }
      if (dim==3)
      {
        out[0][0] = { 0,-1,-1};
        out[0][1] = { 1, 0, 0};
        out[0][2] = { 1, 0, 0};

        out[1][0] = { 0, 1,  0};
        out[1][1] = {-1, 0, -1};
        out[1][2] = { 0, 1,  0};

        out[2][0] = { 0, -1, 0};
        out[2][1] = { 1,  0, 0};
        out[2][2] = { 0,  0, 0};

        out[3][0] = { 0,  0, 1};
        out[3][1] = { 0,  0, 1};
        out[3][2] = {-1, -1, 0};

        out[4][0] = { 0, 0, -1};
        out[4][1] = { 0, 0,  0};
        out[4][2] = { 1, 0,  0};

        out[5][0] = { 0, 0,  0};
        out[5][1] = { 0, 0, -1};
        out[5][2] = { 0, 1,  0};
      }

      for (std::size_t i=0; i<out.size(); i++)
        out[i] *= edgeOrientation_[i];

    }

    /** \brief Evaluate partial derivatives of all shape functions at a given point
     *
     * \param[in] order The partial derivative to be computed, as a multi-index
     * \param[in] in  The evaluation point
     * \param[out] out Jacobians of all shape functions at that point
     */
    void partial(const std::array<unsigned int, dim>& order,
                 const typename Traits::DomainType& in,
                 std::vector<typename Traits::RangeType>& out) const
    {
      auto totalOrder = std::accumulate(order.begin(), order.end(), 0);
      if (totalOrder == 0) {
        evaluateFunction(in, out);
      } else if (totalOrder == 1) {
        auto const direction = std::distance(order.begin(), std::find(order.begin(), order.end(), 1));
        out.resize(size());

        if (dim==2)
        {
          if (direction==0)
          {
            out[0] = {0,  1};
            out[1] = {0, -1};
            out[2] = {0,  1};
          }
          else
          {
            out[0] = {-1, 0};
            out[1] = { 1, 0};
            out[2] = {-1, 0};
          }
        }

        if (dim==3)
        {
          switch (direction)
          {
            case 0:
            out[0] = { 0, 1, 1};
            out[1] = { 0,-1, 0};
            out[2] = { 0, 1, 0};
            out[3] = { 0, 0,-1};
            out[4] = { 0, 0, 1};
            out[5] = { 0, 0, 0};
            break;

            case 1:
            out[0] = {-1, 0, 0};
            out[1] = { 1, 0, 1};
            out[2] = {-1, 0, 0};
            out[3] = { 0, 0,-1};
            out[4] = { 0, 0, 0};
            out[5] = { 0, 0, 1};
            break;

            case 2:
            out[0] = {-1, 0, 0};
            out[1] = { 0,-1, 0};
            out[2] = { 0, 0, 0};
            out[3] = { 1, 1, 0};
            out[4] = {-1, 0, 0};
            out[5] = { 0,-1, 0};
            break;
          }
        }

        for (std::size_t i=0; i<out.size(); i++)
          out[i] *= edgeOrientation_[i];

      } else {
        out.resize(size());
        for (std::size_t i = 0; i < size(); ++i)
          for (std::size_t j = 0; j < dim; ++j)
            out[i][j] = 0;
      }

    }

    //! \brief Polynomial order of the shape functions
    unsigned int order() const
    {
      return k;
    }

  private:

    // Orientations of the simplex edges
    std::array<R,numberOfEdges> edgeOrientation_;
  };


  /** \brief Assignment of Nedelec degrees of freedom to subentities of the reference simplex
   * \tparam dim Dimension of the domain
   * \tparam k Order of the Nedelec finite element
   */
  template <int dim, int k>
  class Nedelec1stKindSimplexLocalCoefficients
  {
  public:
    //! \brief Default constructor
    Nedelec1stKindSimplexLocalCoefficients ()
    : localKey_(size())
    {
      static_assert(k==1, "Only first-order Nédélec local coefficients are implemented.");
      // Assign all degrees of freedom to edges
      // TODO: This is correct only for first-order Nédélec elements
      for (std::size_t i=0; i<size(); i++)
        localKey_[i] = LocalKey(i,dim-1,0);
    }

    //! Number of degrees of freedom
    std::size_t size() const
    {
      static_assert(dim==2 || dim==3, "Nédélec shape functions are implemented only for 2d and 3d simplices.");
      return (dim==2) ? k * (k+2)
                      : k * (k+2) * (k+3) / 2;
    }

    /** \brief Get assignment of i-th degree of freedom to a reference simplex subentity
     */
    const LocalKey& localKey (std::size_t i) const
    {
      return localKey_[i];
    }

  private:
    std::vector<LocalKey> localKey_;
  };

  /** \brief Project a given function into the Nedelec space
   *
   * \tparam LB The LocalBasis that spans the space
   */
  template<class LB>
  class Nedelec1stKindSimplexLocalInterpolation
  {
    static constexpr auto dim = LB::Traits::dimDomain;
    static constexpr auto size = LB::size();

    // Number of edges of the reference simplex
    constexpr static std::size_t numberOfEdges = dim*(dim+1)/2;

  public:

    //! \brief Constructor with given set of edge orientations
    Nedelec1stKindSimplexLocalInterpolation (std::bitset<numberOfEdges> s = 0)
    {
      auto refElement = Dune::referenceElement<double,dim>(GeometryTypes::simplex(dim));

      for (std::size_t i=0; i<numberOfEdges; i++)
        m_[i] = refElement.position(i,dim-1);

      for (std::size_t i=0; i<numberOfEdges; i++)
      {
        auto vertexIterator = refElement.subEntities(i,dim-1,dim).begin();
        auto v0 = *vertexIterator;
        auto v1 = *(++vertexIterator);
        // By default, edges point from the vertex with the smaller index
        // to the vertex with the larger index.
        if (v0>v1)
          std::swap(v0,v1);
        edge_[i] = refElement.position(v1,dim) - refElement.position(v0,dim);
        edge_[i] *= (s[i]) ? -1.0 : 1.0;
      }
    }

    /** \brief Project a given function into the Nedelec space
     *
     * \param f The function to interpolate, must implement RangeField operator(DomainType)
     * \param[out] out The coefficients of the projection
     */
    template<typename F, typename C>
    void interpolate (const F& ff, std::vector<C>& out) const
    {
      out.resize(size);
      auto&& f = Impl::makeFunctionWithCallOperator<typename LB::Traits::DomainType>(ff);

      for (std::size_t i=0; i<size; i++)
      {
        auto y = f(m_[i]);
        out[i] = 0.0;
        for (int j=0; j<dim; j++)
          out[i] += y[j]*edge_[i][j];
      }
    }

  private:
    // Edge midpoints of the reference simplex
    std::array<typename LB::Traits::DomainType, numberOfEdges> m_;
    // Edges of the reference simplex
    std::array<typename LB::Traits::DomainType, numberOfEdges> edge_;
  };

}


  /**
   * \brief Nédélec elements of the first kind for simplex elements
   *
   * These elements have been described in :
   *
   *   J.C. Nédélec, "Mixed finite elements in R^3", Numer.Math., 35(3):315-341,1980.
   *   DOI: http://dx.doi.org/10.1007/BF01396415
   *
   * The order count starts at '1'.  This is the counting used, e.g., by Nédélec himself,
   * and by Kirby, Logg, Rognes, Terrel, "Common and unusual finite elements",
   * https://doi.org/10.1007/978-3-642-23099-8_3
   *
   * \note These shape functions are implemented for the reference simplices only!
   *   The transformation to other simplices has to be done by the user.
   *   One way of doing that is using the implementation of the covariant
   *   Piola transform in globalvaluedlocalfinitelement.hh in dune-functions.
   *
   * \ingroup Nedelec
   *
   * \tparam D Number type used for domain coordinates
   * \tparam R Number type use for shape function values
   * \tparam dim Dimension of the domain
   * \tparam k Order of the Nedelec finite element (lowest is 1)
   *
   */
  template<class D, class R, int dim, int k>
  class Nedelec1stKindSimplexLocalFiniteElement
  {
  public:
    using Traits = LocalFiniteElementTraits<Impl::Nedelec1stKindSimplexLocalBasis<D,R,dim,k>,
                                            Impl::Nedelec1stKindSimplexLocalCoefficients<dim,k>,
                                            Impl::Nedelec1stKindSimplexLocalInterpolation<Impl::Nedelec1stKindSimplexLocalBasis<D,R,dim,k> > >;

    static_assert(dim==2 || dim==3, "Nedelec elements are only implemented for 2d and 3d elements.");
    static_assert(k==1,   "Nedelec elements of the first kind are currently only implemented for order k==1.");

    /** \brief Default constructor
     */
    Nedelec1stKindSimplexLocalFiniteElement() = default;

    /**
     * \brief Constructor with explicitly given edge orientations
     *
     * \param s Edge orientation indicator
     */
    Nedelec1stKindSimplexLocalFiniteElement (std::bitset<dim*(dim+1)/2> s) :
      basis_(s),
      interpolation_(s)
    {}

    const typename Traits::LocalBasisType& localBasis () const
    {
      return basis_;
    }

    const typename Traits::LocalCoefficientsType& localCoefficients () const
    {
      return coefficients_;
    }

    const typename Traits::LocalInterpolationType& localInterpolation () const
    {
      return interpolation_;
    }

    static constexpr unsigned int size ()
    {
      return Traits::LocalBasisType::size();
    }

    static constexpr GeometryType type ()
    {
      return GeometryTypes::simplex(dim);
    }

  private:
    typename Traits::LocalBasisType basis_;
    typename Traits::LocalCoefficientsType coefficients_;
    typename Traits::LocalInterpolationType interpolation_;
  };

}

#endif
