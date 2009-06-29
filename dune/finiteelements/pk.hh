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
  class PkLocalFiniteElement;


  /** \todo Please doc me !
   */
  template<class D, class R>
  class PkLocalFiniteElement<D, R, 1, 1>
    : public P11DLocalFiniteElement<D, R>
  {};

  /** \todo Please doc me !
   */
  template<class D, class R, int k>
  class PkLocalFiniteElement<D, R, 2, k>
    : public Pk2DLocalFiniteElement<D, R, k>
  {};

  /** \todo Please doc me !
   */
  template<class D, class R, int k>
  class PkLocalFiniteElement<D, R, 3, k>
    : public Pk3DLocalFiniteElement<D, R, k>
  {};

}

#endif
