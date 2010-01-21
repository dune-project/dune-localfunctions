// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_REFINED_P0_LOCALBASIS_HH
#define DUNE_REFINED_P0_LOCALBASIS_HH


#include <dune/common/fmatrix.hh>

#include <dune/localfunctions/common/localbasis.hh>

namespace Dune
{
  template<class D, class R, int dim>
  class RefinedP0LocalBasis
    : public RefinedSimplexLocalBasis<D,dim>
  {
  public:
    RefinedP0LocalBasis()
    {
      DUNE_THROW(Dune::NotImplemented,"RefinedP0LocalBasis not implemented for dim > 3.");
    }
  };

  /**@ingroup LocalBasisImplementation
     \brief Uniformly refined constant shape functions on the triangle.

     This shape function set mimicks the P0 shape functions that you would get on
     a uniformly refined grid.  Hence these shape functions are only piecewise
     constant!

     Shape functions like these are necessary for hierarchical error estimators
     for certain nonlinear problems.

     The functions are associated with the simplices having the following centers:

     f_0 ~ (1/6, 1/6)
     f_1 ~ (4/6, 1/6)
     f_2 ~ (1/6, 4/6)
     f_3 ~ (2/6, 2/6)

     \tparam D Type to represent the field in the domain.
     \tparam R Type to represent the field in the range.

     \nosubgrouping
   */
  template<class D, class R>
  class RefinedP0LocalBasis<D,R,2>
    : public RefinedSimplexLocalBasis<D,2>
  {
  public:
    //! \brief export type traits for function signature
    typedef LocalBasisTraits<D,2,Dune::FieldVector<D,2>,R,1,Dune::FieldVector<R,1>, Dune::FieldMatrix<R,1,2> > Traits;

    //! \brief number of shape functions
    unsigned int size () const
    {
      return 4;
    }

    //! \brief Evaluate all shape functions
    inline void evaluateFunction (const typename Traits::DomainType& in,
                                  std::vector<typename Traits::RangeType>& out) const
    {
      out.resize(4);

      switch (getSubElement(in)) {
      case 0 :
        out[0] = 1;
        out[1] = 0;
        out[2] = 0;
        out[3] = 0;
        break;
      case 1 :
        out[0] = 0;
        out[1] = 1;
        out[2] = 0;
        out[3] = 0;
        break;
      case 2 :
        out[0] = 0;
        out[1] = 0;
        out[2] = 1;
        out[3] = 0;
        break;
      case 3 :
        out[0] = 0;
        out[1] = 0;
        out[2] = 0;
        out[3] = 1;
      }
    }

    inline void
    evaluateJacobian (const typename Traits::DomainType& in,         // position
                      std::vector<typename Traits::JacobianType>& out) const      // return value
    {
      out.resize(4);
      out[0][0] = 0;
      out[1][0] = 0;
      out[2][0] = 0;
      out[3][0] = 0;
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
