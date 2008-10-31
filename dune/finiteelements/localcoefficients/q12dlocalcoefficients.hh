// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_Q12DLOCALLAYOUT_HH
#define DUNE_Q12DLOCALLAYOUT_HH

#include <iostream>
#include <vector>

#include "locallayout.hh"

/**@ingroup LocalLayoutImplementation
   \brief Layout map for Q1 elements in 2d.

   \nosubgrouping
 */
class Q12DLocalLayoutMap : public LocalLayoutMapInterface<Q12DLocalLayoutMap>
{
public:
  //! \brief Standard constructor
  Q12DLocalLayoutMap ()
  {
    q12dlayout.resize(4);
    for (int i=0; i<4; i++)
      q12dlayout[i] = LocalIndex(i,2,0);
  }

  //! \brief Deliver layout for entity
  template<class EntityType>
  void find (const EntityType& e, LocalLayout& layout) const
  {
    layout = q12dlayout;
  }

private:
  LocalLayout q12dlayout;
};

#endif
