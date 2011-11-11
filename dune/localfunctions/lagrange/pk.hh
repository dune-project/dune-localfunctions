// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/* vim: set ai expandtab sw=4 ts=4: */
#ifndef DUNE_PK_LOCALFINITEELEMENT_HH
#define DUNE_PK_LOCALFINITEELEMENT_HH

#include "p0.hh"
#include "p1.hh"
#include "pk1d.hh"
#include "pk2d.hh"
#include "pk3d.hh"

namespace Dune
{

  /** \brief General Lagrange finite element with arbitrary dimension and polynomial order
   *
   * \tparam D type used for domain coordinates
   * \tparam R type used for function values
   * \tparam d dimension of the reference element
   * \tparam k polynomial order
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
#if 0
  /** \brief General Lagrange finite element -- specialization for first-order on a 1d reference element
   *
   * \tparam D type used for domain coordinates
   * \tparam R type used for function values
   */
  template<class D, class R>
  class PkLocalFiniteElement<D, R, 1, 1>
    : public P1LocalFiniteElement<D, R, 1>
  {
  public:
    PkLocalFiniteElement()
    {}

    PkLocalFiniteElement(const unsigned int vertexmap[2])
    {}
  };

  /** \brief General Lagrange finite element -- specialization for zero-th order on a 1d reference element
   *
   * \tparam D type used for domain coordinates
   * \tparam R type used for function values
   * \tparam d dimension of the reference element
   */
  template<class D, class R>
  class PkLocalFiniteElement<D, R, 1, 0>
    : public P0LocalFiniteElement<D, R, 1>
  {
  public:
    PkLocalFiniteElement()
      : P0LocalFiniteElement<D,R,1>(GeometryType(0,1))
    {}

    PkLocalFiniteElement(const unsigned int vertexmap[2])
    {}
  };
#endif
  /** \brief General Lagrange finite element -- specialization for a 2d reference element
   *
   * \tparam D type used for domain coordinates
   * \tparam R type used for function values
   * \tparam k polynomial order
   */
  template<class D, class R, int k>
  class PkLocalFiniteElement<D, R, 1, k>
    : public Pk1DLocalFiniteElement<D, R, k>
  {
  public:
    PkLocalFiniteElement()
    {}

    PkLocalFiniteElement(const unsigned int vertexmap[2]) :
      Pk1DLocalFiniteElement<D, R, k>(vertexmap)
    {}
  };

  /** \brief General Lagrange finite element -- specialization for a 2d reference element
   *
   * \tparam D type used for domain coordinates
   * \tparam R type used for function values
   * \tparam k polynomial order
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

  /** \brief General Lagrange finite element -- specialization for a 3d reference element
   *
   * \tparam D type used for domain coordinates
   * \tparam R type used for function values
   * \tparam k polynomial order
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
