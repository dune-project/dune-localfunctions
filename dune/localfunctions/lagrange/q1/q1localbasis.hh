// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_Q1_LOCALBASIS_HH
#define DUNE_Q1_LOCALBASIS_HH

#include <numeric>

#include <dune/common/fmatrix.hh>

#include <dune/localfunctions/common/localbasis.hh>

namespace Dune
{
  /**@ingroup LocalBasisImplementation
         \brief Lagrange shape functions of order 1 on the reference cube.

         Also known as \f$Q^1\f$.

         \tparam D Type to represent the field in the domain.
         \tparam R Type to represent the field in the range.
     \tparam dim Dimension of the cube

         \nosubgrouping
   */
  template<class D, class R, int dim>
  class Q1LocalBasis
  {
  public:
    typedef LocalBasisTraits<D,dim,Dune::FieldVector<D,dim>,R,1,Dune::FieldVector<R,1>,
        Dune::FieldMatrix<R,1,dim> > Traits;

    //! \brief number of shape functions
    unsigned int size () const
    {
      return 1<<dim;
    }

    //! \brief Evaluate all shape functions
    inline void evaluateFunction (const typename Traits::DomainType& in,
                                  std::vector<typename Traits::RangeType>& out) const
    {
      out.resize(size());

      for (size_t i=0; i<size(); i++) {

        out[i] = 1;

        for (int j=0; j<dim; j++)
          // if j-th bit of i is set multiply with in[j], else with 1-in[j]
          out[i] *= (i & (1<<j)) ? in[j] :  1-in[j];

      }

    }

    //! \brief Evaluate Jacobian of all shape functions
    inline void
    evaluateJacobian (const typename Traits::DomainType& in,         // position
                      std::vector<typename Traits::JacobianType>& out) const      // return value
    {
      out.resize(size());

      // Loop over all shape functions
      for (size_t i=0; i<size(); i++) {

        // Loop over all coordinate directions
        for (int j=0; j<dim; j++) {

          // Initialize: the overall expression is a product
          // if j-th bit of i is set to -1, else 1
          out[i][0][j] = (i & (1<<j)) ? 1 : -1;

          for (int k=0; k<dim; k++) {

            if (j!=k)
              // if k-th bit of i is set multiply with in[j], else with 1-in[j]
              out[i][0][j] *= (i & (1<<k)) ? in[k] :  1-in[k];

          }

        }

      }

    }

    /** \brief Evaluate partial derivatives of any order of all shape functions
     * \param order Order of the partial derivatives, in the classic multi-index notation
     * \param in Position where to evaluate the derivatives
     * \param[out] out Return value: the desired partial derivatives
     */
    void partial(const std::array<unsigned int,dim>& order,
                 const typename Traits::DomainType& in,
                 std::vector<typename Traits::RangeType>& out) const
    {
      auto totalOrder = std::accumulate(order.begin(), order.end(), 0);

      if (totalOrder == 0) {
        evaluateFunction(in, out);
      }
      else if (totalOrder == 1) {
        out.resize(size());

        auto direction = std::distance(order.begin(), std::find(order.begin(), order.end(), 1));
        if (direction >= dim) {
          DUNE_THROW(RangeError, "Direction of partial derivative not found!");
        }

        // Loop over all shape functions
        for (std::size_t i = 0; i < size(); ++i) {

          // Initialize: the overall expression is a product
          // if j-th bit of i is set to -1, else 1
          out[i] = (i & (1<<direction)) ? 1 : -1;

          for (int k = 0; k < dim; ++k) {
            if (direction != k)
              // if k-th bit of i is set multiply with in[j], else with 1-in[j]
              out[i] *= (i & (1<<k)) ? in[k] :  1-in[k];
          }

        }
      }
      else
      {
        DUNE_THROW(NotImplemented, "To be implemented!");
      }
    }

    //! \brief Polynomial order of the shape functions
    unsigned int order () const
    {
      return 1;
    }
  };
}
#endif
