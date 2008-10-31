// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_P12DLOCALLAYOUT_HH
#define DUNE_P12DLOCALLAYOUT_HH

#include <iostream>
#include <vector>

#include "locallayout.hh"

/**@ingroup LocalLayoutImplementation
   \brief Layout map for P1 elements in 2d.

   \nosubgrouping
 */
class P12DLocalLayoutMap : public LocalLayoutMapInterface<P12DLocalLayoutMap>
{
public:
  //! \brief Standard constructor
  P12DLocalLayoutMap ()
  {
    p12dlayout.resize(3);
    for (int i=0; i<3; i++)
      p12dlayout[i] = LocalIndex(i,2,0);
  }

  //! \brief Deliver layout for entity
  template<class EntityType>
  void find (const EntityType& e, LocalLayout& layout) const
  {
    layout = p12dlayout;
  }

private:
  LocalLayout p12dlayout;
};

#endif
