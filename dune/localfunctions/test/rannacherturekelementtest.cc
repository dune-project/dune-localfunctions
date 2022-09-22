// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <iostream>

#include <dune/localfunctions/rannacherturek/rannacherturek.hh>

#include <dune/localfunctions/test/test-localfe.hh>

int main(int argc, char** argv)
{
  bool success = true;

  Dune::RannacherTurekLocalFiniteElement<double,double,2> rannacher_turek2dfem;
  TEST_FE(rannacher_turek2dfem);

  Dune::RannacherTurekLocalFiniteElement<double,double,3> rannacher_turek3dfem;
  TEST_FE(rannacher_turek3dfem);

  return success ? 0 : 1;
}
