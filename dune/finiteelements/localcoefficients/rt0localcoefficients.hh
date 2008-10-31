// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_RT0TRIANGLELOCALLAYOUT_HH
#define DUNE_RT0TRIANGLELOCALLAYOUT_HH

#include <iostream>
#include <vector>

#include "locallayout.hh"

/**@ingroup LocalLayoutImplementation
   \brief Layout map for RT0 elements on triangles.

   \nosubgrouping
 */
class RT0TriangleLocalLayoutMap : public LocalLayoutMapInterface<RT0TriangleLocalLayoutMap>
{
public:

  //! \brief Standard constructor
  RT0TriangleLocalLayoutMap ()
  {
    rt0layout.resize(3);
    for (int i=0; i<3; i++)
      rt0layout[i] = LocalIndex(i,1,0);
  }

  //! \brief Deliver layout for entity
  template<class EntityType>
  void find (const EntityType& e, LocalLayout& layout) const
  {
    layout = rt0layout;
  }

private:
  LocalLayout rt0layout;
};

#endif
