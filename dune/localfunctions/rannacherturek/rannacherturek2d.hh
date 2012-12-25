// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_RANNACHER_TUREK_2D_LOCALFINITEELEMENT_HH
#define DUNE_RANNACHER_TUREK_2D_LOCALFINITEELEMENT_HH

#ifndef DISABLE_RANNACHERTUREK2D_DEPRECATION_WARNING
#warning dune/localfunctions/rannacherturek/rannacherturek2d.hh\
  is deprecated, please use dune/localfunctions/rannacherturek/rannacherturek.hh
#endif

#include "rannacherturek.hh"

namespace Dune
{
  /**
   * \deprecated This class is deprecated and will be removed after Dune 2.3.
   *             Use RannacherTurekLocalFiniteElement< D, R, 2 > instead.
   */
  template< class D, class R >
  class
  DUNE_DEPRECATED_MSG("Use RannacherTurekLocalFiniteElement< D, R, 2 > instead")
  RannacherTurek2DLocalFiniteElement
    : public RannacherTurekLocalFiniteElement< D, R, 2 >
  {};

} // namespace Dune

#endif // #ifndef DUNE_RANNACHER_TUREK_2D_LOCALFINITEELEMENT_HH
