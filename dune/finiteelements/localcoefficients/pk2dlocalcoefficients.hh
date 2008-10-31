// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PK2DLOCALLAYOUT_HH
#define DUNE_PK2DLOCALLAYOUT_HH

#include <iostream>
#include <vector>

#include "locallayout.hh"

/**@ingroup LocalLayoutImplementation
   \brief Layout map for Pk elements in 2d.

   \nosubgrouping
 */
template<typename GV, unsigned int k>
class Pk2DLocalLayoutMap : public LocalLayoutMapInterface<Pk2DLocalLayoutMap<GV,k> >
{
  enum {N = (k+1)*(k+2)/2};
  typedef typename GV::IndexSet IndexSet;

public:

  //! \brief Standard constructor
  Pk2DLocalLayoutMap (const GV& gv_) : gv(gv_), is(gv_.indexSet())
  {
    pk2dlayout.resize(N);
    int n=0;
    int c=0;
    for (int j=0; j<=k; j++)
      for (int i=0; i<=k-j; i++)
      {
        if (n==0)
        {
          pk2dlayout[n++] = LocalIndex(0,2,0);
          continue;
        }
        if (i==k && j==0)
        {
          pk2dlayout[n++] = LocalIndex(1,2,0);
          continue;
        }
        if (i==0 && j==k)
        {
          pk2dlayout[n++] = LocalIndex(2,2,0);
          continue;
        }
        if (j==0)
        {
          pk2dlayout[n++] = LocalIndex(2,1,i-1);
          continue;
        }
        if (i==0)
        {
          pk2dlayout[n++] = LocalIndex(1,1,j-1);
          continue;
        }
        if (i+j==k)
        {
          pk2dlayout[n++] = LocalIndex(0,1,j-1);
          continue;
        }
        pk2dlayout[n++] = LocalIndex(0,0,c);
        c++;
      }
  }

  //! \brief Deliver layout for entity
  template<class EntityType>
  void find (const EntityType& e, LocalLayout& layout) const
  {
    layout = pk2dlayout;
    if (k<3) return;

    unsigned int n0,n1,n2;
    n0 = is.template subIndex<2>(e,0);
    n1 = is.template subIndex<2>(e,1);
    n2 = is.template subIndex<2>(e,2);
    bool flip[3];
    flip[0] = n1>n2;
    flip[1] = n0>n2;
    flip[2] = n0>n1;
    for (int i=0; i<layout.size(); i++)
      if (layout[i].codim()==1 && flip[layout[i].subentity()])
        layout[i].index(k-2-layout[i].index());

    //  for (int i=0; i<layout.size(); i++)
    //    std::cout << "n=" << i << " s=" << layout[i].subentity()
    //                          << " c=" << layout[i].codim()
    //                          << " i=" << layout[i].index();
  }

private:
  const GV& gv;
  const IndexSet& is;
  LocalLayout pk2dlayout;
};


template<typename GV>
class Pk2DLocalLayoutMap<GV,0> : public LocalLayoutMapInterface<Pk2DLocalLayoutMap<GV,0> >
{
public:
  Pk2DLocalLayoutMap (const GV& gv)
  {
    p0layout.resize(1);
    p0layout[0] = LocalIndex(0,0,0);
  }

  template<class EntityType>
  void find (const EntityType& e, LocalLayout& layout) const
  {
    layout = p0layout;
  }

private:
  LocalLayout p0layout;
};

#endif
