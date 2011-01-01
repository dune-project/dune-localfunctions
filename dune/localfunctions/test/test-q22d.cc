// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <cstddef>
#include <iostream>
#include <ostream>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/geometrytype.hh>

#include <dune/grid/utility/mockgeometry.hh>

#include <dune/localfunctions/lagrange/q22d.hh>

#include "test-fe.hh"

int main(int argc, char** argv) {
  try {
    // tolerance for floating-point comparisons
    static const double eps = 1e-9;
    // stepsize for numerical differentiation
    static const double delta = 1e-5;

    int result = 77;

    {
      std::cout << "== Checking global-valued Q22D elements" << std::endl;

      Dune::GeometryType gt;
      gt.makeQuadrilateral();

      Dune::FieldVector<double, 2> corners[4];
      corners[0][0] = -.5; corners[0][1] = 0;
      corners[1][0] = 0  ; corners[1][1] = -.5;
      corners[2][0] =  .5; corners[2][1] = 0;
      corners[3][0] = 0  ; corners[3][1] =  .5;
      typedef Dune::MockGeometry<double, 2, 2> Geometry;
      Geometry geo(gt, corners);

      Dune::Q22DFiniteElementFactory<Geometry, double> feFactory;
      bool success = testFE(geo, feFactory.make(geo), eps, delta);

      if(success && result != 1)
        result = 0;
      else
        result = 1;
    }

    return result;
  }
  catch (const Dune::Exception& e) {
    std::cerr << e << std::endl;
    throw;
  }
}
