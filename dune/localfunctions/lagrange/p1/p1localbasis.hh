// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_P1_LOCALBASIS_HH
#define DUNE_P1_LOCALBASIS_HH

#include <array>
#include <numeric>

#include <dune/common/deprecated.hh>
#include <dune/common/fmatrix.hh>

#include <dune/localfunctions/common/localbasis.hh>

namespace Dune
{
  /**@ingroup LocalBasisImplementation
         \brief Linear Lagrange shape functions on the simplex.

         Defines the linear shape functions on the simplex.

         \tparam D Type to represent the field in the domain.
         \tparam R Type to represent the field in the range.
     \tparam dim The dimension of the simplex

         \nosubgrouping
   */
  template<class D, class R, int dim>
  class P1LocalBasis
  {
  public:
    //! \brief export type traits for function signature
    typedef LocalBasisTraits<D,dim,Dune::FieldVector<D,dim>,R,1,Dune::FieldVector<R,1>,
        Dune::FieldMatrix<R,1,dim>, 2> Traits;

    //! \brief number of shape functions
    unsigned int size () const
    {
      return dim+1;
    }

    //! \brief Evaluate all shape functions
    inline void evaluateFunction (const typename Traits::DomainType& in,
                                  std::vector<typename Traits::RangeType>& out) const
    {
      out.resize(size());
      out[0] = 1.0;
      for (size_t i=0; i<dim; i++) {
        out[0]  -= in[i];
        out[i+1] = in[i];
      }
    }

    //! \brief Evaluate Jacobian of all shape functions
    inline void
    evaluateJacobian (const typename Traits::DomainType& in,         // position
                      std::vector<typename Traits::JacobianType>& out) const      // return value
    {
      out.resize(size());

      for (int i=0; i<dim; i++)
        out[0][0][i] = -1;

      for (int i=0; i<dim; i++)
        for (int j=0; j<dim; j++)
          out[i+1][0][j] = (i==j);

    }

    /** \brief Evaluate partial derivatives of any order of all shape functions
     * \param order Order of the partial derivatives, in the classic multi-index notation
     * \param in Position where to evaluate the derivatives
     * \param[out] out Return value: the desired partial derivatives
     */
    inline void partial(const std::array<unsigned int,dim>& order,
                        const typename Traits::DomainType& in,
                        std::vector<typename Traits::RangeType>& out) const
    {
      auto totalOrder = std::accumulate(order.begin(), order.end(), 0);

      if (totalOrder==0)
        evaluateFunction(in, out);
      else if (totalOrder==1)
      {
        auto direction = std::find(order.begin(), order.end(), 1);
        out.resize(size());

        out[0] = -1;
        for (int i=0; i<dim; i++)
          out[i+1] = (i==(direction-order.begin()));
      }
      else  // all higher order derivatives are zero
      {
        out.resize(size());

        for (int i=0; i<dim+1; i++)
          out[i] = 0;
      }
    }

    //! \brief Evaluate all shape functions
    template<unsigned int k>
    inline void DUNE_DEPRECATED_MSG("Use method 'partial' instead!")
    evaluate (const typename std::array<int,k>& directions,
                          const typename Traits::DomainType& in,
                          std::vector<typename Traits::RangeType>& out) const
    {
      if (k==0)
        evaluateFunction(in, out);
      else if (k==1)
      {
        out.resize(size());

        out[0] = -1;
        for (int i=0; i<dim; i++)
          out[i+1] = (i==directions[0]);
      }
      else if (k==2)
      {
        out.resize(size());

        for (int i=0; i<dim+1; i++)
          out[i] = 0;
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
