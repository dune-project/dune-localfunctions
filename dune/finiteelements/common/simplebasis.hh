// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_SIMPLEBASIS_HH
#define DUNE_SIMPLEBASIS_HH

#include <vector>

#include <dune/common/static_assert.hh>
#include <dune/common/typetraits.hh>

#include "basis.hh"

namespace Dune
{

  template<typename LB>
  class C0SimpleBasis;
  /**@ingroup LocalBasisInterface
   * \brief BasisTraits for C0SimpleBasis
   *
   * \tparam LB  Type of the local basis.
   * \tparam Geo Type of the geometry.
   *
   * \nosubgrouping
   */
  template<typename LB>
  struct BasisTraits<C1SimpleBasis<LB> >
    : public typename LB::Traits
  {
    //! Limit the diffOrder to 0
    static const unsigned diffOrder = 0;
  };
  /**@ingroup LocalBasisInterface
   * \brief Basis implementation where local and global values are equal
   *
   * It simply takes a localbasis and returns whatever that returns.  This is
   * most appropriate for single component functions.  For multicomponent
   * functions there is typically some transformation needed.  Which
   * transformation to use depends very much on the particular base function.
   *
   * \tparam LB Type of the local basis.
   *
   * \nosubgrouping
   */
  template<typename LB>
  class C0SimpleBasis
    : public C0BasisInterface<C0SimpleBasis<LB> >
  {
    const LB& lb;

  public:
    //! \brief Export type traits
    typedef BasisTraits<C0SimpleBasis<LB> > Traits;

    // construct a C0SimpleBasis from a localbasis
    /**
     * \param [in] lb_ The localbasis to wrap.  A reference to this object is
     *                 stored, that means, the object must live at least as
     *                 long as this simple basis is used for evaluation.
     */
    C0SimpleBasis(const LB& lb_)
      : lb(lb_)
    {}

    //! \copydoc C0BasisInterface::size
    unsigned int size () const
    {
      return lb.size();
    }

    //! \copydoc C0BasisInterface::evaluateFunction
    inline void
    evaluateFunction(const typename Traits::DomainType& in,
                     std::vector<typename Traits::RangeType>& out) const
    {
      lb.evaluateFunction(in,out);
    }

    //! \copydoc C0BasisInterface::order
    unsigned int order () const
    {
      return lb.order();
    }
  };

  template<typename LB, typename Geo>
  class C1SimpleBasis;
  /**@ingroup LocalBasisInterface
   *  \brief BasisTraits for C1SimpleBasis
   * \tparam LB  Type of the local basis.
   * \tparam Geo Type of the geometry.
   *
   * \nosubgrouping
   */
  template<typename LB, typename Geo>
  struct BasisTraits<C1SimpleBasis<LB, Geo> >
    : public typename LB::Traits
  {
    dune_static_assert(Geo::mydimension == dimDomain,
                       "Local dimension of the geometry and "
                       "domain dimension of the local basis must match");
    dune_static_assert(LB::Traits::diffOrder >= 1,
                       "diffOrder of local basis must be at least 1 "
                       "for C1SimpleBasis");

    //! Limit the diffOrder to 1
    static const unsigned diffOrder = 1;
    //! The number of columns of the jacobian is different than for the LB
    typedef FieldMatrix<
        typename LB::Traits::JacobianType::field_type::field_type,
        dimRange, Geo::coorddimension> JacobianType;
  };
  /**@ingroup LocalBasisInterface
   * \brief Basis implementation where local and global values are equal
   *
   * It simply takes a localbasis and returns whatever that returns.  This is
   * most appropriate for single component functions.  For multicomponent
   * functions there is typically some transformation needed.  Which
   * transformation to use depends very much on the particular base function.
   *
   * \tparam LB  Type of the local basis.
   * \tparam Geo Type of the geometry.  Its interface should conform to that
   *             of Dune::Geometry.
   *
   * \nosubgrouping
   */
  template<typename LB, typename Geo>
  class C1SimpleBasis
    : public C1BasisInterface<C1SimpleBasis<LB, Geo> >
  {
    const LB& lb;
    const Geo& geo;

  public:
    //! \brief Export type traits
    typedef BasisTraits<C1SimpleBasis<LB, Geo> > Traits;

    //! construct a C1SimpleBasis from a localbasis
    /**
     * \param [in] lb_ The localbasis to wrap.  A reference to this object is
     *                 stored, that means, the object must live at least as
     *                 long as this simple basis is used for evaluation.
     */
    C1SimpleBasis(const LB& lb_, const Geo& geo_)
      : lb(lb_), geo(geo_)
    {}

    //! \copydoc C1BasisInterface::size
    unsigned int size () const
    {
      return lb.size();
    }

    //! \copydoc C1BasisInterface::evaluateFunction
    inline void
    evaluateFunction(const typename Traits::DomainType& in,
                     std::vector<typename Traits::RangeType>& out) const
    {
      lb.evaluateFunction(in,out);
    }

    //! \copydoc C1BasisInterface::order
    unsigned int order () const
    {
      return lb.order();
    }

    //! \copydoc C1BasisInterface::evaluateJacobian
    /**
     *  This implementation is for single component base functions.  For
     *  multicomponent base functions, it will transform the gradient of each
     *  component independently, as if each component was an independent
     *  scalar valued base function.
     *
     * The Jacobian is taken as a vector of the gradients of each component:
     * \f[
     *    J_f=\begin{pmatrix}
     *          (\nabla f_0)^T\\
     *          (\nabla f_1)^T\\
     *          \vdots
     *        \end{pmatrix}
     * \f]
     * Each gradient is transformed according to
     * \f[
     *    \nabla\phi^i_\alpha|_{\hat x} =
     *       J^{-T}_g|_{\hat x} \hat\nabla\hat\phi^i_\alpha|_{\hat x},
     * \f]
     * where \f$\phi^i_\alpha\f$ is component \f$\alpha\f$ of the \f$i\f$'th
     * shape function, \f$J_g\f$ is the Jacobian of the local-to-global map of
     * the geometry and the hat \f$\hat{}\f$ denotes local coordinates
     * (\f$\hat x\f$) or functions of the local basis (as opposed to this
     * "global" basis).
     */
    inline void
    evaluateJacobian(const typename Traits::DomainType& in,             // position
                     std::vector<typename Traits::JacobianType>& out) const                          // return value
    {
      out.resize(lb.size());
      std::vector<typename LB::Traits::JacobianType> localJ(lb.size());
      lb.evaluateJacobian(in, localjacobian);
      typename Geo::Jacobian geoJinvT = geo.jacobianInverseTransposed(in);
      for(unsigned baseno = 0; baseno < lb.size(); ++baseno)
        for(unsigned compno = 0; compno < Traits::dimRange; ++compno)
          geoJinvT.mv(localJ[baseno][compno], out[baseno][compno]);
    }

  };
}
#endif // DUNE_SIMPLEBASIS_HH
