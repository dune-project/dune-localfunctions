// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_RANNACHER_TUREK_2D_LOCALCOEFFICIENTS_HH
#define DUNE_RANNACHER_TUREK_2D_LOCALCOEFFICIENTS_HH

#warning dune/localfunctions/rannacherturek/rannacherturek2d/rannacherturek2dlocalcoefficients.h\
  is deprecated, please use dune/localfunctions/rannacherturek/rannachertureklocalcoefficients.hh

#include <dune/localfunctions/rannacherturek/rannachertureklocalcoefficients.hh>

namespace Dune
{

  class RannacherTurek2DLocalCoefficients
    : public RannacherTurekLocalCoefficients< 2 >
  {};

} //namespace Dune

#endif // #ifndef DUNE_RANNACHER_TUREK_2D_LOCALCOEFFICIENTS_HH
