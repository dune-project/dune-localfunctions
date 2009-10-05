// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef BASISPRINT
#define BASISPRINT
#include <dune/finiteelements/multiindex.hh>
namespace Dune {
  template <class Basis>
  void basisPrint(Basis &basis)
  {
    const int dimension = Basis::dimension;
    typedef MultiIndex< dimension > Field;

    unsigned int size = basis.size();

    std::cout << "Polynomial representation of the basis functions:" << std::endl;
    std::cout << "Number of base functions:  " << size << std::endl;
    std::vector< Dune::FieldVector<Field,dimension> > y( size );
    FieldVector< Field, dimension > x;
    for( int i = 0; i < dimension; ++i )
      x[ i ].set( i, 1 );
    basis.evaluate( x, y );
    std::cout << y << std::endl;
  }
};
#endif // BASISPRINT
