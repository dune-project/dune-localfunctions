// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DGORTHONORMALBASIS_HH
#define DUNE_DGORTHONORMALBASIS_HH

#include <dune/finiteelements/orthonormalbasis/orthonormalbasis.hh>
#include <dune/finiteelements/dglocalcoefficients.hh>

namespace Dune
{

  template< int dim, class SF, class CF = typename ComputeField< SF, 512 >::Type >
  class DGOrthonormalBasisProvider
    : public BasisProvider< ONBasisCreator< dim, SF, CF > >
  {
    typedef DGOrthonormalBasisProvider< dim, SF, CF > This;
    typedef BasisProvider< ONBasisCreator< dim, SF, CF > > Base;

  public:
    static const unsigned int dimension = dim;

    typedef typename Base::Key Key;

    typedef DGLocalCoefficients LocalCoefficients;

    template< class Topology >
    static const LocalCoefficients &localCoefficients ( const Key &key )
    {
      const unsigned int numCoefficients = Base::basis( Topology::id, key ).size();
      return *(new LocalCoefficients( numCoefficients ));
    }

    static void release ( const LocalCoefficients &localCoefficients )
    {
      delete &localCoefficients;
    }
  };

}

#endif // #ifndef DUNE_DGORTHONORMALBASIS_HH
