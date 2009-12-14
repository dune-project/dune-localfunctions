// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_REFINED_P0_LOCALBASIS_HH
#define DUNE_REFINED_P0_LOCALBASIS_HH


#include "../common/localbasis.hh"

namespace Dune
{
  template<class D, class R, int dim>
  class RefinedP0LocalBasis :
    public C1LocalBasisInterface<
        C1LocalBasisTraits<D,dim,Dune::FieldVector<D,dim>,R,1,Dune::FieldVector<R,1>, Dune::FieldVector<Dune::FieldVector<R,dim>,1> >
#ifndef DUNE_VIRTUAL_SHAPEFUNCTIONS
        ,RefinedP0LocalBasis<D,R,dim>
#endif
        >
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
  class RefinedP0LocalBasis<D,R,2> :
    public C1LocalBasisInterface<
        C1LocalBasisTraits<D,2,Dune::FieldVector<D,2>,R,1,Dune::FieldVector<R,1>, Dune::FieldVector<Dune::FieldVector<R,2>,1> >
        // We derive from C1 to allow using this with the virtual interface.
        // Otherwise the convariant return type of the localBasis() method
        // in the FE does not match the base class.
        // Notice that this is not even C0.
#ifndef DUNE_VIRTUAL_SHAPEFUNCTIONS
        ,RefinedP0LocalBasis<D,R,2>
#endif
        >
  {
  public:
    //! \brief export type traits for function signature
    typedef C1LocalBasisTraits<D,2,Dune::FieldVector<D,2>,R,1,Dune::FieldVector<R,1>, Dune::FieldVector<Dune::FieldVector<R,2>,1> > Traits;

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

  private:
    /** \brief Get local coordinates in the subtriangle.
     *
     * The triangles are ordered according to
     *
     * |\
     * |2\
     * |--\
     * |\3|\
     * |0\|1\
     * ------
     *
     * \param[in] global Coordinates in the reference triangle
     * \returns Number of the subtriangles containing in
     */
    static int getSubElement(const typename Traits::DomainType& global)
    {
      if (global[0] + global[1] <= 0.5)
        return 0;
      else if (global[0] >= 0.5)
        return 1;
      else if (global[1] >= 0.5)
        return 2;
      return 3;
    }
  };

}
#endif
