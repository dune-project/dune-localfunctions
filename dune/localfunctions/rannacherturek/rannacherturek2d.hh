// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_RANNACHER_TUREK_2D_LOCALFINITEELEMENT_HH
#define DUNE_RANNACHER_TUREK_2D_LOCALFINITEELEMENT_HH

#warning dune/localfunctions/rannacherturek/rannacherturek2d.hh\
  is deprecated, please use\
  dune/localfunctions/rannacherturek/rannacherturek.hh

#include "rannacherturek.hh"

namespace Dune
{

  template< class D, class R >
  class RannacherTurek2DLocalFiniteElement
    : public RannacherTurekLocalFiniteElement< D, R, 2 >
  {};

} // namespace Dune

#endif // #ifndef DUNE_RANNACHER_TUREK_2D_LOCALFINITEELEMENT_HH
