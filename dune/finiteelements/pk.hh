// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/* vim: set ai expandtab sw=4 ts=4: */
#ifndef DUNE_PK_LOCALFINITEELEMENT_HH
#define DUNE_PK_LOCALFINITEELEMENT_HH

#include "p11d.hh"
#include "pk2d.hh"
#include "pk3d.hh"


namespace Dune
{

  /** \todo Please doc me !
   */
  template<class D, class R, int d, int k>
  class PkLocalFiniteElement
  {
  public:
    PkLocalFiniteElement()
    {}

    /** Constructor for variants with permuted vertices.

        \param vertexmap The permutation of the vertices.  This
        can for instance be generated from the global indices of
        the vertices by reducing those to the integers 0...k+1
     */
    PkLocalFiniteElement(const unsigned int vertexmap[k+1])
    {}
  };

  /** \todo Please doc me !
   */
  template<class D, class R>
  class PkLocalFiniteElement<D, R, 1, 1>
    : public P11DLocalFiniteElement<D, R>
  {
  public:
    PkLocalFiniteElement()
    {}

    PkLocalFiniteElement(const unsigned int vertexmap[2])
    {}
  };

  /** \todo Please doc me !
   */
  template<class D, class R, int k>
  class PkLocalFiniteElement<D, R, 2, k>
    : public Pk2DLocalFiniteElement<D, R, k>
  {
  public:
    PkLocalFiniteElement()
    {}

    PkLocalFiniteElement(const unsigned int vertexmap[3]) :
      Pk2DLocalFiniteElement<D, R, k>(vertexmap)
    {}
  };

  /** \todo Please doc me !
   */
  template<class D, class R, int k>
  class PkLocalFiniteElement<D, R, 3, k>
    : public Pk3DLocalFiniteElement<D, R, k>
  {
  public:
    PkLocalFiniteElement()
    {}

    PkLocalFiniteElement(const unsigned int vertexmap[4]) :
      Pk3DLocalFiniteElement<D, R, k>(vertexmap)
    {}
  };

}

#endif
