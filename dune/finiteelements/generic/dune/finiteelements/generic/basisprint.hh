// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef BASISPRINT
#define BASISPRINT
#include <dune/finiteelements/generic/multiindex.hh>
namespace Dune {
  template <int deriv,class Basis>
  void basisPrint(std::ostream &out, Basis &basis)
  {
    const int dimension = Basis::dimension;
    typedef MultiIndex< dimension > Field;

    unsigned int size = basis.size();

    out << "% Number of base functions:  " << size << std::endl;
    out << "% Derivative order: " << deriv << std::endl;
    std::vector< Dune::FieldVector<Field,Basis::dimRange> > y( size );
    FieldVector< Field, dimension > x;
    for( int i = 0; i < dimension; ++i )
      x[ i ].set( i, 1 );
    basis.template evaluate<deriv>( x, y );
    out << y << std::endl;
  }
};
#endif // BASISPRINT
