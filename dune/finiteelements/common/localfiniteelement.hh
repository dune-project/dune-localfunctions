// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LOCALFINITEELEMENT_HH
#define DUNE_LOCALFINITEELEMENT_HH

#include <iostream>
#include <vector>

#include <dune/common/geometrytype.hh>

#include <dune/finiteelements/common/localbasis.hh>
#include <dune/finiteelements/common/localcoefficients.hh>
#include <dune/finiteelements/common/localinterpolation.hh>

namespace Dune {

  //! traits helper struct
  template<class LB, class LC, class LI>
  struct LocalFiniteElementTraits
  {
    /** \todo Please doc me !
     */
    typedef LB LocalBasisType;

    /** \todo Please doc me !
     */
    typedef LC LocalCoefficientsType;

    /** \todo Please doc me !
     */
    typedef LI LocalInterpolationType;
  };

  /** \brief interface for a finite element
      \tparam T The controlling traits
   */
#if DUNE_VIRTUAL_SHAPEFUNCTIONS
  template <class DT, class RT, int dim>
  class LocalFiniteElementInterface
  {
  public:

    class Traits
    {

      typedef Dune::FieldVector<RT,1> RangeType;

    };

    typedef Dune::C1LocalBasisInterface<Dune::C1LocalBasisTraits<DT,
            dim,
            Dune::FieldVector<DT,dim>,
            RT,
            1,
            Dune::FieldVector<RT,1>,
            Dune::FieldVector<Dune::FieldVector<RT,dim>,1> > > LocalBasisType;

    /** \todo Please doc me !
     */
    const virtual LocalBasisType& localBasis () const = 0;

    /** \todo Please doc me !
     */
    const virtual LocalCoefficientsInterface& localCoefficients () const = 0;

    /** \todo Please doc me !
     */
    const virtual LocalInterpolationInterface& localInterpolation () const = 0;

    /** \todo Please doc me !
     */
    virtual GeometryType type () const = 0;

  };

#else
  template<class T, class Imp>
  class LocalFiniteElementInterface
  {
  public:

    /** \todo Please doc me !
     */
    typedef T Traits;

    /** \todo Please doc me !
     */
    const typename T::LocalBasisType& localBasis () const
    {
      return asImp().localBasis();
    }

    /** \todo Please doc me !
     */
    const typename T::LocalCoefficientsType& localCoefficients () const
    {
      return asImp().localCoefficients();
    }

    /** \todo Please doc me !
     */
    const typename T::LocalInterpolationType& localInterpolation () const
    {
      return asImp().localInterpolation();
    }

    /** \todo Please doc me !
     */
    GeometryType type () const
    {
      return asImp().type();
    }

  private:
    Imp& asImp () {return static_cast<Imp &> (*this);}
    const Imp& asImp () const {return static_cast<const Imp &>(*this);}
  };
#endif

}

#endif
