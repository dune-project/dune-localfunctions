// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <cstddef>
#include <iostream>

#include <dune/localfunctions/common/virtualinterface.hh>

#include <dune/localfunctions/p1.hh>

/** \file
    \brief Test the dynamically polymorphic shape function interface
 */

using namespace Dune;

const int dim = 2;

typedef double D;
typedef double R;

//typedef C0LocalBasisTraits<D,dim,FieldVector<D,dim>,R,1,FieldVector<R,1> > C0Traits;
typedef C1LocalBasisTraits<D,2,FieldVector<D,2>,R,1,FieldVector<R,1>,
                           FieldVector<FieldVector<R,2>,1> > C1Traits;

struct MagicLocalFiniteElementTraits
{
  typedef C1LocalBasisVirtualInterface<C1Traits> LocalBasisType;

};


int main (int argc, char *argv[]) try
{

  const Dune::LocalFiniteElementVirtualImp<MagicLocalFiniteElementTraits,
      Dune::P1LocalFiniteElement<D, R, dim> > virtualLocalSimplexFE_;

  return 0;

}
catch (Exception e) {

  std::cout << e << std::endl;
  return 1;
}
