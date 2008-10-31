// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_Q22DLOCALLAYOUT_HH
#define DUNE_Q22DLOCALLAYOUT_HH

#include <iostream>
#include <vector>

#include "locallayout.hh"

/**@ingroup LocalLayoutImplementation
   \brief Layout map for Q2 elements in 2d.

   \nosubgrouping
 */
class Q22DLocalLayoutMap : public LocalLayoutMapInterface<Q22DLocalLayoutMap>
{
public:
  //! \brief Standard constructor
  Q22DLocalLayoutMap ()
  {
    q22dlayout.resize(9);
    for (int i=0; i<9; i++)
      q22dlayout[i] = LocalIndex(i%4,2-i/4,0);
  }

  //! \brief Deliver layout for entity
  template<class EntityType>
  void find (const EntityType& e, LocalLayout& layout) const
  {
    layout = q22dlayout;
  }

private:
  LocalLayout q22dlayout;
};

#endif
