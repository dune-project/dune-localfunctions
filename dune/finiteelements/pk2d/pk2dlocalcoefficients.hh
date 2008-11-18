// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PK2DLOCALCOEFFICIENTS_HH
#define DUNE_PK2DLOCALCOEFFICIENTS_HH

#include <iostream>
#include <vector>

#include "../common/localcoefficients.hh"

namespace Dune
{

  /**@ingroup LocalLayoutImplementation
         \brief Layout map for P0 elements

         \nosubgrouping
   */
  template<unsigned int k>
  class Pk2DLocalCoefficients : public LocalCoefficientsInterface<Pk2DLocalCoefficients<k> >
  {
    enum {N = (k+1)*(k+2)/2};

  public:
    //! \brief Standard constructor
    Pk2DLocalCoefficients () : li(N)
    {
      fill_default();
    }

    //! constructor for eight variants with order on edges flipped
    Pk2DLocalCoefficients (int variant) : li(N)
    {
      fill_default();
      bool flip[3];
      if (variant & 1) flip[0]=true;else flip[0]=false;
      if (variant & 2) flip[1]=true;else flip[1]=false;
      if (variant & 4) flip[2]=true;else flip[2]=false;
      for (int i=0; i<N; i++)
        if (li[i].codim()==1 && flip[li[i].subentity()])
          li[i].index(k-2-li[i].index());
    }

    //! number of coefficients
    int size () const
    {
      return N;
    }

    //! get i'th index
    const LocalIndex& localIndex (int i) const
    {
      return li[i];
    }

  private:
    std::vector<LocalIndex> li;

    void fill_default ()
    {
      if (k==0)
      {
        li[0] = LocalIndex(0,0,0);
        return;
      }
      int n=0;
      int c=0;
      for (unsigned int j=0; j<=k; j++)
        for (unsigned int i=0; i<=k-j; i++)
        {
          if (i==0 && j==0)
          {
            li[n++] = LocalIndex(0,2,0);
            continue;
          }
          if (i==k && j==0)
          {
            li[n++] = LocalIndex(1,2,0);
            continue;
          }
          if (i==0 && j==k)
          {
            li[n++] = LocalIndex(2,2,0);
            continue;
          }
          if (j==0)
          {
            li[n++] = LocalIndex(2,1,i-1);
            continue;
          }
          if (i==0)
          {
            li[n++] = LocalIndex(1,1,j-1);
            continue;
          }
          if (i+j==k)
          {
            li[n++] = LocalIndex(0,1,j-1);
            continue;
          }
          li[n++] = LocalIndex(0,0,c++);
        }
    }
  };

}

#endif
