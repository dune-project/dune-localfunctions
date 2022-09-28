// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALKEY_HH
#define DUNE_LOCALKEY_HH

#include <array>
#include <cstddef>
#include <ostream>

namespace Dune
{
  /**@ingroup LocalLayoutInterface
         \brief Describe position of one degree of freedom

         A LocalKey associates a degree of freedom with an index
         of a local basis function.

         \nosubgrouping
   */
  class LocalKey
  {
  public:

    /** \brief Enumerate 'special values' for the codimension method */
    enum {
      /** \brief Codimension returned by LocalKey::codim() for degrees of freedom attached to an intersection

         The standard interface of dune-localfunctions assumes that degrees of freedom are attached to subentities
         of an element.  This subentities can be described by a codimension and a subentity number.
         However some elements, like the mimetic finite elements, attach their degrees of freedom to intersections.
         While intersections do have a codimension, namely 1, having the method codim() return 1 in this case
         would be ambiguous.  Hence 'intersectionCodim' is returned instead.
       */
      intersectionCodim=666
    };

    //! \brief Standard constructor for uninitialized local index
    LocalKey ()
    {}

    /** \brief Initialize all components
        \param s Local number of the associated subentity
        \param c Codimension of the associated subentity
        \param i Index in the set of all functions associated to this subentity
     */
    LocalKey (unsigned int s, unsigned int c, unsigned int i)
    {
      values_[0] = s;
      values_[1] = c;
      values_[2] = i;
    }

    //! \brief Return number of associated subentity
    inline unsigned int subEntity () const
    {
      return values_[0];
    }

    //! \brief Return codim of associated entity
    inline unsigned int codim () const
    {
      return values_[1];
    }

    //! \brief Return offset within subentity
    inline unsigned int index () const
    {
      return values_[2];
    }

    //! \brief Set index component
    void index (unsigned int i)
    {
      values_[2] = i;
    }

    /** \brief Less-than operator so we can use this class as a key type in stl containers */
    bool operator< (const LocalKey& other) const
    {
      return values_ < other.values_;
    }

    /** \brief Write LocalKey object to output stream */
    friend std::ostream& operator<< (std::ostream& s, const LocalKey& localKey)
    {
      return s << "[ subEntity: " << localKey.subEntity()
             << ", codim: " << localKey.codim()
             << ", index: " << localKey.index() << " ]";
    }

  private:

    // We use an array to store the values in order to be able to use the array::operator< implementation
    std::array<unsigned int,3> values_;

  };

}
#endif
