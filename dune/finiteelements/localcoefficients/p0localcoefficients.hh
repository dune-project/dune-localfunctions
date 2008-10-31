// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_P0LOCALLAYOUT_HH
#define DUNE_P0LOCALLAYOUT_HH

#include <iostream>
#include <vector>

#include "locallayout.hh"

/**@ingroup LocalLayoutImplementation
   \brief Layout map for P0 elements

   \nosubgrouping
 */
class P0LocalLayoutMap : public LocalLayoutMapInterface<P0LocalLayoutMap>
{
public:
  //! \brief Standard constructor
  P0LocalLayoutMap ()
  {
    p0layout.resize(1);
    p0layout[0] = LocalIndex(0,0,0);
  }

  //! \brief Deliver layout for entity
  template<class EntityType>
  void find (const EntityType& e, LocalLayout& layout) const
  {
    layout = p0layout;
  }

private:
  LocalLayout p0layout;
};

#endif
