// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LOCALFINITEELEMENT_HH
#define DUNE_LOCALFINITEELEMENT_HH

#include <iostream>
#include <vector>

#include <dune/common/interfaces.hh>
#include <dune/common/geometrytype.hh>

#include <dune/localfunctions/common/localbasis.hh>
#include <dune/localfunctions/common/localcoefficients.hh>
#include <dune/localfunctions/common/localinterpolation.hh>

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

#if DUNE_VIRTUAL_SHAPEFUNCTIONS
  /** \brief interface for a finite element
      \tparam DT The type used for domain coordinates
      \tparam RT The type used for range coordinates
      \tparam dim Domain dimension
   */
  template <class DT, class RT, int dim>
  class LocalFiniteElementInterface :
    public Cloneable
  {
  public:
    /** \brief Export various types.  The weird nested structure
        is there to mimick the static inheritance case.
     */
    struct Traits
    {
      typedef typename Dune::C1LocalBasisInterface<Dune::C1LocalBasisTraits<
              DT,
              dim,
              Dune::FieldVector<DT,dim>,
              RT,
              1,
              Dune::FieldVector<RT,1>,
              Dune::FieldVector<Dune::FieldVector<RT,dim>,1> > > LocalBasisType;

      typedef Dune::LocalCoefficientsInterface LocalCoefficientsType;
      typedef typename Dune::LocalInterpolationInterface<typename LocalBasisType::Traits::DomainType, typename LocalBasisType::Traits::RangeType> LocalInterpolationType;

    };

    /** \todo Please doc me !
     */
    const virtual typename Traits::LocalBasisType& localBasis () const = 0;

    /** \todo Please doc me !
     */
    const virtual typename Traits::LocalCoefficientsType& localCoefficients () const = 0;

    /** \todo Please doc me !
     */
    const virtual typename Traits::LocalInterpolationType& localInterpolation () const = 0;

    /** \todo Please doc me !
     */
    virtual GeometryType type () const = 0;

    virtual LocalFiniteElementInterface* clone () const = 0;

  };

#else
  /** \brief interface for a finite element
      \tparam T The controlling traits
      \tparam Imp The deriving implementation class
   */
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

    //      Imp* clone () const
    //      {
    //        return asImp().clone();
    //      }

  private:
    Imp& asImp () {return static_cast<Imp &> (*this);}
    const Imp& asImp () const {return static_cast<const Imp &>(*this);}
  };
#endif

}

#endif
