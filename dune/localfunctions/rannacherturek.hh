// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/** \file
    \brief Convenience header that includes all available Rannacher-Turek LocalFiniteElements
 */


#include <dune/localfunctions/rannacherturek/rannacherturek.hh>

#define DISABLE_RANNACHERTUREK2D_DEPRECATION_WARNING
#include <dune/localfunctions/rannacherturek/rannacherturek2d.hh>
#undef DISABLE_RANNACHERTUREK2D_DEPRECATION_WARNING
