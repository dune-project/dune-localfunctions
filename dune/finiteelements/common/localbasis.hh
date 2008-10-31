// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LOCALBASIS_HH
#define DUNE_LOCALBASIS_HH

#include <iostream>
#include <vector>

#include <dune/common/static_assert.hh>
#include <dune/common/fixedarray.hh>
#include <dune/common/fvector.hh>

/**@ingroup LocalBasisInterface
   \brief Type traits for LocalBasisInterface

   A shape function is a function
   \f[ \hat\phi : \mbox{IR}^n \to \mbox{IR}^m. \f]
   This traits class holds information how the signature of this
   function is represented in C++ types.

   This is just a convenience class for supplying traits to the
   C0LocalBasisInterface.

   Template parameters:

   - <tt>DF</tt>: Type to represent the field in the domain.
   - <tt>n</tt>:  Dimension of the domain.
   - <tt>D</tt>:  Type to represent the domain, allows random access.
   - <tt>RF</tt>: Type to represent the field in the range.
   - <tt>m</tt>:  Dimension of the range.
   - <tt>R</tt>:  Type to represent the range, allows random access.

   \nosubgrouping
 */
template<class DF, int n, class D, class RF, int m, class R>
struct C0LocalBasisTraits
{
  //! \brief Export type for domain field
  typedef DF DomainFieldType;

  //! \brief Enum for domain dimension
  enum {
    //! \brief dimension of the domain
    dimDomain = n
  };

  //! \brief domain type
  typedef D DomainType;

  //! \brief Export type for range field
  typedef RF RangeFieldType;

  //! \brief Enum for range dimension
  enum {
    //! \brief dimension of the range
    dimRange = m
  };

  //! \brief range type
  typedef R RangeType;

  //! \brief Enum for differentiability order
  enum {
    //! \brief number of derivatives supported
    diffOrder=0
  };
};

/**@ingroup LocalBasisInterface
   \brief Interface for shape functions on a specific reference element

   This class represents a set of shape functions defined on one particular
   reference element.

   Template parameters:

   - <tt>T</tt>: Instance of LocalBasisTraits providing type information.
   - <tt>Imp</tt>: Implementation of the interface used via Barton-Nackman

   \nosubgrouping
 */
template<class T, class Imp>
class C0LocalBasisInterface
{
public:
  //! \brief Export type traits
  typedef T Traits;

  //! \brief Number of shape functions
  unsigned int size () const
  {
    return asImp().size();
  }

  /** \brief Evaluate all basis function at given position

      Evaluates all shape functions at the given position and returns
      these values in a vector.
   */
  inline void evaluateFunction (const typename Traits::DomainType& in,
                                std::vector<typename Traits::RangeType>& out) const
  {
    asImp().evaluateFunction(in,out);
  }

  /** \brief Represent given function with the shape functions

      Determines a set of coefficients returned in out
      to represent a given function on an element.
      If the function is not contained within the local function
      spaced spanned by the shape functions an approximation
      is returned.

      This function is to be moved to the Layout because it depends on the way
      how the basis functions are used!
   */
  template<typename E, typename F, typename C>
  void interpolate (const E& e, const F& f, std::vector<C>& out) const
  {
    asImp().interpolate(e,f,out);
  }

  /*! \brief Polynomial order of the shape functions

     \todo Gurke!
   */
  unsigned int order () const
  {
    return asImp().order();
  }

private:
  Imp& asImp () {return static_cast<Imp &> (*this);}
  const Imp& asImp () const {return static_cast<const Imp &>(*this);}
};



/**@ingroup LocalBasisInterface
   \brief Type traits for C1LocalBasisInterface

   Extends the traits class LocalBasisTraits for differentiable
   shape functions.

   Template parameters:

   - <tt>DF</tt>: Type to represent the field in the domain.
   - <tt>n</tt>:  Dimension of the domain.
   - <tt>D</tt>:  Type to represent the domain, allows random access.
   - <tt>RF</tt>: Type to represent the field in the range.
   - <tt>m</tt>:  Dimension of the range.
   - <tt>R</tt>:  Type to represent the range, allows random access.
   - <tt>J</tt>:  Type to represent the Jacobian, allows random access.

   \nosubgrouping
 */
template<class DF, int n, class D, class RF, int m, class R, class J>
struct C1LocalBasisTraits : public C0LocalBasisTraits<DF,n,D,RF,m,R>
{
  /** \brief Type to represent derivative

      When \f$ \hat\phi : \mbox{IR}^n \to \mbox{IR}^m \f$ then JacobianType
      is an 2D-array of m x n components where entry J[i][j] contains
      the derivative  \f$\partial_j \hat\phi_i \f$.
   */
  typedef J JacobianType;

  //! \brief Enum for differentiability order
  enum {
    //! \brief number of derivatives supported
    diffOrder=1
  };
};



/**@ingroup LocalBasisInterface
   \brief Interface for differentiable shape functions on a specific reference element

   This class represents a set of differentiable shape functions defined on one particular
   reference element.

   Template parameters:

   - <tt>T</tt>: Instance of C1LocalBasisTraits providing type information.
   - <tt>Imp</tt>: Implementation of the interface used via Barton-Nackman

   \nosubgrouping
 */
template<class T, class Imp>
class C1LocalBasisInterface : public C0LocalBasisInterface<T,Imp>
{
public:
  //! \brief Export type traits
  typedef T Traits;

  /** \brief Evaluate jacobian of all shape functions at given position.

      out[k][i][j] is \f$\partial_i \hat\phi_j^k \f$, when \f$\hat\phi^k \f$ is the
      k'th shape function.
   */
  inline void
  evaluateJacobian(const typename Traits::DomainType& in,         // position
                   std::vector<typename Traits::JacobianType>& out) const      // return value
  {
    asImp().evaluateJacobian(in,out);
  }

private:
  Imp& asImp () {return static_cast<Imp &> (*this);}
  const Imp& asImp () const {return static_cast<const Imp &>(*this);}
};



template<class DF, int n, class D, class RF, int m, class R, class J, int dorder>
struct CkLocalBasisTraits : public C1LocalBasisTraits<DF,n,D,RF,m,R,J>
{
  //! \brief Enum for differentiability order
  enum {
    //! \brief number of derivatives supported
    diffOrder=dorder
  };
};

// hinzufügen:
template<class T, class Imp>
class CkLocalBasisInterface : public C1LocalBasisInterface<T,Imp>
{
public:
  //! \brief Export type traits
  typedef T Traits;

  //! return given derivative of all components
  template<int k>
  inline void evaluate (const Dune::array<int,k>& directions,
                        const typename Traits::DomainType& in,
                        std::vector<typename Traits::RangeType>& out) const
  {
    asImp().evaluate(directions,in,out);
  }

private:
  Imp& asImp () {return static_cast<Imp &> (*this);}
  const Imp& asImp () const {return static_cast<const Imp &>(*this);}
};

#endif
