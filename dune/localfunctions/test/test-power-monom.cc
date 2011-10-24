// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <cstddef>
#include <iostream>
#include <ostream>

#include <dune/common/exceptions.hh>
#include <dune/common/forloop.hh>

#include <dune/localfunctions/meta/power.hh>
#include <dune/localfunctions/monom.hh>

#include "geometries.hh"
#include "test-fe.hh"

template<int dimD>
struct DimD {
  template<int dimR>
  struct DimR {
    template<int p>
    struct Order {
      // tolerance for floating-point comparisons
      // be a little bit more lenient here than usual, the momom local basis
      // is known to become more and more unstable with increasing order
      static constexpr double eps = 1e-8;
      // stepsize for numerical differentiation
      static constexpr double delta = 1e-5;

      static void apply(int &result) {
        std::cout << "== Checking global-valued Power elements (with "
                  << "dimR=" << dimR << ") wrapping Monom elements (with "
                  << "dimD=" << dimD << ", p=" << p << ")" << std::endl;

        typedef TestGeometries<double, dimD> TestGeos;
        static const TestGeos testGeos;

        typedef typename TestGeos::Geometry Geometry;
        typedef Dune::MonomFiniteElementFactory<Geometry, double, p>
        BackendFEFactory;
        BackendFEFactory backendFEFactory;

        typedef typename BackendFEFactory::FiniteElement BackendFE;
        Dune::PowerFiniteElementFactory<BackendFE, dimR> feFactory;

        for(std::size_t i = 0; i < testGeos.size(); ++i) {
          const Geometry &geo = testGeos[i];

          std::cout << "=== GeometryType " << geo.type() << std::endl;

          bool success =
            testFE(geo, feFactory.make(backendFEFactory.make(geo)), eps,
                   delta);

          if(success && result != 1)
            result = 0;
          else
            result = 1;
        }
      }
    };

    static void apply(int &result)
    { Dune::ForLoop<Order, 0, 3>::apply(result); }
  };

  static void apply(int &result)
  { Dune::ForLoop<DimR, 0, 3>::apply(result); }
};

int main(int argc, char** argv) {
  try {
    int result = 77;

    Dune::ForLoop<DimD, 1, 3>::apply(result);

    return result;
  }
  catch (const Dune::Exception& e) {
    std::cerr << e << std::endl;
    throw;
  }
}
