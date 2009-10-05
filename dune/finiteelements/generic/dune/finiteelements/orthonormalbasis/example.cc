// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// http://www.alglib.net/eigen/symmetric/symmevd.php
// Uebersetzen:
// g++ -O3 -o example example.cc -Ilibs -lgmp -lmpfr -I$HOME/MORGOTH/dune/dune-fem -I$HOME/MORGOTH/dune/dune-common
// g++ -O3 -o example example.cc -Ilibs libmpfr.a libgmp.a
const unsigned int Precision = 1024;
#include "config.h"
#include "orthonormalcompute.hh"
using namespace Dune;
template <int dim,int ord>
void printCoeff(std::ostream& out,
                mat_t& res) {
  int N = res.gethighbound(1);
  for (int i=1; i<=N; ++i) {
    out << "Polynomial : " << i << std::endl;
    for (int j=1; j<=i; j++) {
      if (fabs(res(j,i).toDouble())<1e-20)
        out << 0 << "\t\t" << std::flush;
      else
        out << amp::ampf<128>(res(j,i)).toDec() << "\t\t" << std::flush;
    }
    for (int j=i+1; j<=N; j++) {
      assert(fabs(res(j,i).toDouble())<1e-10);
    }
    out << std::endl;
  }
}
/******************************************/
template <int dim,int ord>
void cube(CalcCoeffs<dim,ord>& calc) {
  // cube
  MultiIndex<dim> geo;
  std::stringstream name;
  name << "cube" << dim << "d.coef";
  for (int i=0; i<dim; ++i)
    geo.set(i,2);
  calc.compute(geo);
  std::ofstream cube(name.str().c_str());
  printCoeff<dim,ord>(cube,calc.res);
}
template <int dim,int ord>
void simplex(CalcCoeffs<dim,ord>& calc) {
  // simplex
  MultiIndex<dim> geo;
  std::stringstream name;
  name << "simplex" << dim << "d.coef";
  for (int i=0; i<dim; ++i)
    geo.set(i,1);
  calc.compute(geo);
  std::ofstream cube(name.str().c_str());
  printCoeff<dim,ord>(cube,calc.res);
}
template <int dim,int nr,int ord>
struct AllCubeSimplex {
  static void calc() {
    {
      mat_t res;
      cube<dim,ord>(res);
      simplex<dim,ord>(res);
    }
    AllCubeSimplex<dim+1,nr-1,ord>::calc();
  }
};
template <int dim,int ord>
struct AllCubeSimplex<dim,0,ord> {
  static void calc() {
    {
      mat_t res;
      cube<dim,ord>(res);
      simplex<dim,ord>(res);
    }
  }
};
/******************************************/
/******************************************/
/******************************************/
int main(int argc, char ** argv, char ** env) {
  int use_method = atoi(argv[1]);
  const int ord = 10;
  std::stringstream name;
  std::string geotype;
  {
    std::stringstream classname;
    CalcCoeffs<1,ord> calc(use_method);
    cube<1,ord>(calc);
    std::ofstream cub("cube1d.coef");
    printCoeff<1,ord>(cub,calc.res);
  }
  { // 2d:
    CalcCoeffs<2,ord> calc(use_method);
    {
      cube<2,ord>(calc);
      std::ofstream cub("cube2d.coef");
      printCoeff<2,ord>(cub,calc.res);
    }
    {
      simplex<2,ord>(calc);
      std::ofstream sim("simplex2d.coef");
      printCoeff<2,ord>(sim,calc.res);
    }
  }
  { // 3d:
    CalcCoeffs<3,ord> calc(use_method);
    {
      cube<3,ord>(calc);
      std::ofstream cub("cube3d.coef");
      printCoeff<3,ord>(cub,calc.res);
    }
    {
      simplex<3,ord>(calc);
      std::ofstream sim("simplex3d.coef");
      printCoeff<3,ord>(sim,calc.res);
    }
    MultiIndex<3> geo;
    { // prisma
      geo.set(0,1);
      geo.set(1,1);
      geo.set(2,2);
      calc.compute(geo);
      std::ofstream prism("prism3d.coef");
      printCoeff<3,ord>(prism,calc.res);
    }
    { // pyramid
      geo.set(0,2);
      geo.set(1,2);
      geo.set(2,1);
      calc.compute(geo);
      std::ofstream prism("pyramid3d.coef");
      printCoeff<3,ord>(prism,calc.res);
    }
  }
}
