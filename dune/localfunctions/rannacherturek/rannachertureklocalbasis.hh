// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_RANNACHER_TUREK_LOCALBASIS_HH
#define DUNE_RANNACHER_TUREK_LOCALBASIS_HH

#include "rannacherturek2d/rannacherturek2dlocalbasis.hh"
#include "rannacherturek3d/rannacherturek3dlocalbasis.hh"

namespace Dune
{

  /**@ingroup LocalBasisImplementation
     \brief Rannacher-Turek shape functions

     \tparam D type to represent the field in the domain.
     \tparam R type to represent the field in the range.
     \tparam d domain dimension

     \nosubgrouping
   */
  template< class D, class R, unsigned int d >
  struct RannacherTurekLocalBasis;

  template< class D, class R >
  struct RannacherTurekLocalBasis< D, R, 2 >
    : public RannacherTurek2DLocalBasis< D, R >
  {};

  template< class D, class R >
  struct RannacherTurekLocalBasis< D, R, 3 >
    : public RannacherTurek3DLocalBasis< D, R >
  {};

} // namespace Dune

#endif // #ifndef DUNE_RANNACHER_TUREK_LOCALBASIS_HH
