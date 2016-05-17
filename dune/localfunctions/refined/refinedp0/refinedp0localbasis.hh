// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_REFINED_P0_LOCALBASIS_HH
#define DUNE_REFINED_P0_LOCALBASIS_HH

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <dune/localfunctions/common/localbasis.hh>
#include <dune/localfunctions/refined/common/refinedsimplexlocalbasis.hh>

namespace Dune
{

  /**@ingroup LocalBasisImplementation
     \brief Uniformly refined constant shape functions on a unit simplex in R^dim

     This shape function set mimicks the P0 shape functions that you would get on
     a uniformly refined grid.  Hence these shape functions are only piecewise
     constant!

     Shape functions like these are necessary for hierarchical error estimators
     for certain nonlinear problems.

     The functions are associated with the subelements as defined in RefinedSimplexLocalBasis

     \tparam D Type to represent the field in the domain.
     \tparam R Type to represent the field in the range.
     \tparam dim Dimension of domain space

     \nosubgrouping
   */
  template<class D, class R, int dim>
  class RefinedP0LocalBasis
    : public RefinedSimplexLocalBasis<D,dim>
  {
    // 2 to the k-th power
    enum {N = 1<<dim};
  public:
    //! \brief export type traits for function signature
    typedef LocalBasisTraits<D,dim,Dune::FieldVector<D,dim>,R,1,
      Dune::FieldVector<R,1>, Dune::FieldMatrix<R,1,dim>, DUNE_MAX_DIFF_ORDER > Traits;

    //! \brief number of shape functions
    constexpr std::size_t size () const
    {
      return N;
    }

    //! \brief Evaluate all shape functions
    inline void evaluateFunction (const typename Traits::DomainType& in,
                                  std::vector<typename Traits::RangeType>& out) const
    {
      int subElement = this->getSubElement(in);
      out.resize(N);
      for(int i=0; i<N; ++i)
        out[i] = (i==subElement) ? 1 : 0;
    }

    inline void
    evaluateJacobian (const typename Traits::DomainType& /*in*/,         // position
                      std::vector<typename Traits::JacobianType>& out) const      // return value
    {
      out.resize(N);
      for(int i=0; i<N; ++i)
        out[i][0] = 0;
    }

    //! \brief Evaluate partial derivatives of all shape functions
    inline void partial (const std::array<unsigned int, dim>& order,
                         const typename Traits::DomainType& in,         // position
                         std::vector<typename Traits::RangeType>& out) const      // return value
    {
      auto totalOrder = std::accumulate(order.begin(), order.end(), 0);
      if (totalOrder == 0) {
        evaluateFunction(in, out);
      } else {
        out.resize(size());
        for (std::size_t i = 0; i < size(); ++i)
          out[i] = 0;
      }
    }

    //! \brief Evaluate partial derivatives of all shape functions
    template <std::size_t dOrder>
    inline void evaluate (const std::array<int, dOrder>& /*directions*/,
                          const typename Traits::DomainType& in,         // position
                          std::vector<typename Traits::RangeType>& out) const      // return value
    {
      if (dOrder == 0) {
        evaluateFunction(in, out);
      } else {
        out.resize(size());
        for (std::size_t i = 0; i < size(); ++i)
          out[i] = 0;
      }
    }

    /** \brief Polynomial order of the shape functions
     *
     * Doesn't really apply: these shape functions are only piecewise constant
     */
    unsigned int order () const
    {
      return 0;
    }

  };

}
#endif
