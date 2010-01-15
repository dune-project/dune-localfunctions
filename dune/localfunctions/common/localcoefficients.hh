// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LOCALCOEFFICIENTS_HH
#define DUNE_LOCALCOEFFICIENTS_HH

#include <cstddef>

#include <dune/common/tuples.hh>

namespace Dune
{
  // the number of the beast ...
  enum {
    //! codim that indicates degree of freedom in intersection
    intersectionCodim=666
  };

  /**@ingroup LocalLayoutInterface
         \brief Describe position of one degree of freedom

         A LocalKey associates a degree of freedom with an index
         of a local basis function.

         \nosubgrouping
   */
  class LocalKey
  // The data members of this class are implemented using a tuple base class,
  // because that way they can be used as keys in stl containers.
    : public Dune::tuple<unsigned int, unsigned int, unsigned int>
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

}
#endif
