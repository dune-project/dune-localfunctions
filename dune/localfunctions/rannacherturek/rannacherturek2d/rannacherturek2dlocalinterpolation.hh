// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_RANNACHER_TUREK_2D_LOCALINTERPOLATION_HH
#define DUNE_RANNACHER_TUREK_2D_LOCALINTERPOLATION_HH

#warning dune/localfunctions/rannacherturek/rannacherturek2d/rannacherturek2dlocalinterpolation.hh\
  is deprecated, please use dune/localfunctions/rannacherturek/rannachertureklocalinterpolation.hh

#include <dune/localfunctions/rannacherturek/rannachertureklocalinterpolation.hh>

namespace Dune
{

  template< class LB >
  class RannacherTurek2DLocalInterpolation
    : public RannacherTurekLocalInterpolation< typename LB::Traits::DomainFieldType,
          typename LB::Traits::RangeFieldType,
          LB::Traits::dimRange
          >
  {};

} // namespace Dune

#endif // #ifndef DUNE_RANNACHER_TUREK_2D_LOCALINTERPOLATION_HH
