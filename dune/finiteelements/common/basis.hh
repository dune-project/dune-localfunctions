// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_BASIS_HH
#define DUNE_BASIS_HH

#include <vector>

#include <dune/common/static_assert.hh>

namespace Dune
{

  //! Traits for basis implementations
  /**
   *  Each basis should either have its own specialization of this class or
   *  provide a traits class conforming to this interface in some other way.
   *
   *  \note Defining the traits inside the basis class is impossible if the
   *  basis is derived from C0BasisInterface or C1BasisInterface since those
   *  require the traits class as an (implicit or explicit) template
   *  parameter.
   */
  template<typename Imp>
  struct BasisTraits {
    dune_static_assert(false,
                       "If you get this error you the compiler tried to "
                       "instantiate the non-specialized version of "
                       "BasisTraits, which should never happen.  Maybe you "
                       "forgot to provide a specialization for your basis?");
    //! Type used for single coordinate components
    typedef Imp::Traits::DomainFieldType DomainFieldType;
    //! Dimension of the domain, number of components per coordinate
    static const unsigned dimDomain = Imp::Traits::dimDomain;
    //! Type used for complete coordinates
    typedef Imp::Traits::DomainType DomainType;

    //! Type used for one component of the function value
    typedef Imp::Traits::RangeFieldType RangeFieldType;
    //! Dimension or number of components of the function value
    static const unsigned dimRange = Imp::Traits::dimRange;
    //! Type used for the complete function value
    typedef Imp::Traits::RangeType RangeType;

    //! How many times this function may be differentiated
    /**
     * This should be 0 if the basis provides only the interface of
     * C0BasisInterface and 1 if it provides the interface of
     * C1BasisInterface.
     */
    static const unsigned diffOrder Imp::Traits::diffOrder;
    //! The type of the jacobian
    /**
     *  \note This typedef is not required for a basis that implements the
     *        C0BasisInterface only
     */
    typedef Imp::Traits::JacobianType JacobianType;
  };

  /**@ingroup LocalBasisInterface
         \brief Interface for shape functions on a specific reference element

         This class represents a set of shape functions defined on one particular
     reference element.  It returns global values.

     \tparam Imp Implementation of the interface used via CRTP
     \tparam T   Instance of LocalBasisTraits providing type information.

         \nosubgrouping
   */
  template<typename Imp, typename T = BasisTraits<Imp> >
  class C0BasisInterface
  {
  public:
    //! \brief Export type traits
    typedef T Traits;

    //! \brief Number of shape functions
#if DUNE_VIRTUAL_SHAPEFUNCTIONS
    virtual unsigned int size () const = 0;
#else
    unsigned int size () const
    {
      return asImp().size();
    }
#endif

    /** \brief Evaluate all basis function at given position
     *
     *  Evaluates all shape functions at the given (local) position and
     *  returns the (global) values in a vector.
     *
     *  \param [in]  in  Where to evaluate in local (reference element)
     *                   coordinates.
     *  \param [out] out The resulting global values, one per shape function.
     */
#if DUNE_VIRTUAL_SHAPEFUNCTIONS
    virtual void evaluateFunction(const typename Traits::DomainType& in,
                                  std::vector<typename Traits::RangeType>& out) const = 0;
#else
    inline void evaluateFunction(const typename Traits::DomainType& in,
                                 std::vector<typename Traits::RangeType>& out) const
    {
      asImp().evaluateFunction(in,out);
    }
#endif

    /** \brief Evaluate basis with coefficients at given position
     *
     *  Evaluates all shape functions at the given (local) position and
     *  return the (global) value of the sum weighted with some coefficients.
     *
     *  \note This is a default implementation which works by perfoming the
     *        weighted sum of the values returned by evaluateFunction().
     *
     *  \tparam C The type of the coefficients
     *
     *  \param [in] in     Where to evaluate in local (reference element)
     *                     coordinates.
     *  \param [in] coeffs The coefficients.
     *  \return     The resulting global value.
     */
    template<typename C>
    inline typename Traits::RangeType
    evaluateCoeffs(const typename Traits::DomainType& in,
                   const std::vector<C>& coeffs) const
    {
      std::vector<typename Traits::RangeType> basevalues;
      typename Traits::RangeType out = 0;
      asImp().evaluateFunction(in, basevalues);
      assert(coeffs.size() == basevalues.size());
      for(unsigned i = 0; i < coeffs.size(); ++i)
        out.axpy(coeffs[i], basevalues[i]);
      return out;
    }

    /*! \brief Polynomial order of the shape functions

       \todo Gurke!
     */
#if DUNE_VIRTUAL_SHAPEFUNCTIONS
    virtual unsigned int order () const = 0;
#else
    unsigned int order () const
    {
      return asImp().order();
    }
#endif

#ifndef DUNE_VIRTUAL_SHAPEFUNCTIONS
  private:
    Imp& asImp () {return static_cast<Imp &> (*this);}
    const Imp& asImp () const {return static_cast<const Imp &>(*this);}
#endif
  };



  /**@ingroup LocalBasisInterface
   * \brief Interface for differentiable shape functions on a specific
   *        reference element
   *
   * This class represents a set of differentiable shape functions defined on
   * one particular reference element.  This interface returns global values.
   *
   * \tparam Imp Implementation of the interface used via CRTP.
   * \tparam T   Instance of C1LocalBasisTraits providing type information.
   *
   * \nosubgrouping
   */
  template<typename Imp, typename T = BasisTraits<Imp> >
  class C1BasisInterface
    : public C0BasisInterface<Imp,T>
  {
  public:
    //! \brief Export type traits
    typedef T Traits;

    /** \brief Evaluate jacobian of all shape functions at given position.
     *
     *  out[i][j][k] is \f$\partial_k \hat\phi_j^i \f$, when \f$\hat\phi^i \f$
     *  is the i'th shape function.  Note that this are the derivatives of the
     *  global values by the global coordinates, evaluated at local
     *  coordinates.
     *
     *  \param [in]  in  Where to evaluate in local coordinates.
     *  \param [out] out The result, one jacobian per base function.
     */
#if DUNE_VIRTUAL_SHAPEFUNCTIONS
    virtual void
    evaluateJacobian(const typename Traits::DomainType& in,             // position
                     std::vector<typename Traits::JacobianType>& out) const = 0;
#else
    inline void
    evaluateJacobian(const typename Traits::DomainType& in,             // position
                     std::vector<typename Traits::JacobianType>& out) const                          // return value
    {
      asImp().evaluateJacobian(in,out);
    }
#endif

    /** \brief Evaluate jacobian with given coefficients at given position.
     *
     *  out[j][k] is \f$\partial_k \hat f_j \f$, where \f$f\f$ is the function
     *  given by this basis and the coefficients.  Note that this are the
     *  derivatives of the global values by the global coordinates, evaluated
     *  at local coordinates.
     *
     *  \note This is a default implementation which works by perfoming the
     *        weighted sum of the values returned by evaluateJacobian().
     *
     *  \tparam C The type of the coefficients
     *
     *  \param [in] in     Where to evaluate in local (reference element)
     *                     coordinates.
     *  \param [in] coeffs The coefficients.
     *  \return            The resulting global jacobian.
     */
    template<typename C>
    inline void
    evaluateJacobianCoeffs(const typename Traits::DomainType& in,         // position
                           const std::vector<C>& coeffs) const      // return value
    {
      std::vector<typename Traits::JacobianType> basevalues;
      typename Traits::JacobianType out = 0;
      asImp().evaluateFunction(in, basevalues);
      assert(coeffs.size() == basevalues.size());
      for(unsigned i = 0; i < coeffs.size(); ++i)
        out.axpy(coeffs[i], basevalues[i]);
      return out;
    }

#ifndef DUNE_VIRTUAL_SHAPEFUNCTIONS
  private:
    Imp& asImp () {return static_cast<Imp &> (*this);}
    const Imp& asImp () const {return static_cast<const Imp &>(*this);}
#endif
  };
}
#endif // DUNE_BASIS_HH
