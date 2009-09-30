// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_BASIS_HH
#define DUNE_BASIS_HH

#include <vector>

#include <dune/common/array.hh>
#include <dune/common/fvector.hh>

namespace Dune
{

  /**@ingroup LocalBasisInterface
         \brief Interface for shape functions on a specific reference element

         This class represents a set of shape functions defined on one particular
     reference element.  It returns global values.

         \tparam T     Instance of LocalBasisTraits providing type information.
     \tparam Imp   Implementation of the interface used via CRTP

         \nosubgrouping
   */
  template<class T
#ifndef DUNE_VIRTUAL_SHAPEFUNCTIONS
      , class Imp
#endif
      >
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
     *  \tparam C The type of the coefficients
     *
     *  \param [in]  in     Where to evaluate in local (reference element)
     *                      coordinates.
     *  \param [in]  coeffs The coefficients.
     *  \param [out] out    The resulting global value.
     */
    template<typename C>
    inline void evaluateCoeffs(const typename Traits::DomainType& in,
                               const std::vector<C>& coeffs,
                               typename Traits::RangeType& out) const
    {
      asImp().evaluateCoeffs(in,coeffs,out);
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
   * \tparam T   Instance of C1LocalBasisTraits providing type information.
   * \tparam Imp Implementation of the interface used via CRTP.
   *
   * \nosubgrouping
   */
#if DUNE_VIRTUAL_SHAPEFUNCTIONS
  template<class T>
  class C1BasisInterface
    : public C0BasisInterface<T>
#else
  template<class T, class Imp>
  class C1BasisInterface
    : public C0BasisInterface<T,Imp>
#endif
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
     *  \tparam C The type of the coefficients
     *
     *  \param [in]  in     Where to evaluate in local (reference element)
     *                      coordinates.
     *  \param [in]  coeffs The coefficients.
     *  \param [out] out    The resulting global jacobian.
     */
    template<typename C>
    inline void
    evaluateJacobianCoeffs(const typename Traits::DomainType& in,         // position
                           const std::vector<C>& coeffs,
                           typename Traits::JacobianType& out) const      // return value
    {
      asImp().evaluateJacobianCoeffs(in,coeffs,out);
    }

#ifndef DUNE_VIRTUAL_SHAPEFUNCTIONS
  private:
    Imp& asImp () {return static_cast<Imp &> (*this);}
    const Imp& asImp () const {return static_cast<const Imp &>(*this);}
#endif
  };
}
#endif // DUNE_BASIS_HH
