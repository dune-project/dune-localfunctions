// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LOCALCOEFFICIENTS_HH
#define DUNE_LOCALCOEFFICIENTS_HH

#include <iostream>
#include <vector>

#include <dune/common/tuples.hh>

namespace Dune
{
  /**@ingroup LocalLayoutInterface
         \brief Describe position of one degree of freedom

         A LocalKey associates a degree of freedom with an index
         of a local basis function.

         \nosubgrouping
   */
  class LocalKey : public Dune::tuple<unsigned int, unsigned int, unsigned int>
  {
  public:
    //! \brief Standard constructor for uninitialized local index
    LocalKey ()
    {}

    /** \brief Initialize all components
        \param s Local number of the associated subentity
        \param c Codimension of the associated subentity
        \param i Index in the set of all functions associated to this subentity
     */
    LocalKey (unsigned int s, unsigned int c, unsigned int i)
      : Dune::tuple<unsigned int, unsigned int, unsigned int>(s,c,i)
    {}

    //! \brief Return number of associated subentity
    inline unsigned int subEntity () const
    {
      return Dune::get<0>(*this);
    }

    //! \brief Return codim of associated entity
    inline unsigned int codim () const
    {
      return Dune::get<1>(*this);
    }

    //! \brief Return offset within subentity
    inline unsigned int index () const
    {
      return Dune::get<2>(*this);
    }

    //! \brief Set index component
    void index (unsigned int i)
    {
      Dune::get<2>(*this) = i;
    }
  };

  /**@ingroup LocalLayoutInterface
         \brief Layout description corresponding to a local basis

         To make a (global) finite element function out of a
         local basis on each element certain basis functions in touching
         elements have to be identified.

         This is achieved by mapping the index of a local basis function to a geometric
         (sub-)entity and an offset within that entity.

         \nosubgrouping
   */
#ifndef DUNE_VIRTUAL_SHAPEFUNCTIONS
  template<class Imp>
#endif
  class LocalCoefficientsInterface
  {
  public:
    //! number of coefficients
#if DUNE_VIRTUAL_SHAPEFUNCTIONS
    virtual int size () const = 0;
#else
    int size () const
    {
      return asImp().size();
    }
#endif

    //! get i'th index
#if DUNE_VIRTUAL_SHAPEFUNCTIONS
    const virtual LocalKey& localKey (int i) const = 0;
#else
    const LocalKey& localKey (int i) const
    {
      return asImp().localKey(i);
    }
#endif

#ifndef DUNE_VIRTUAL_SHAPEFUNCTIONS
  private:
    Imp& asImp () {return static_cast<Imp &> (*this);}
    const Imp& asImp () const {return static_cast<const Imp &>(*this);}
#endif
  };

}
#endif
