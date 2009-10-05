// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ORTHONORMALBASIS_HH
#define DUNE_ORTHONORMALBASIS_HH
#include <sstream>
#include "../coeffmatrix.hh"
#include "../monomialbasis.hh"
#include "../multiindex.hh"
#include "orthonormalcompute.hh"
namespace Dune
{
  template <class Topology,class SF,class CF=ComputeField<SF,512> >
  class OrthonormalBasis
  {
    enum {dimension = Topology::dimension};

    typedef typename SF StorageField;
    typedef typename CF ComputationField;

    typedef StandardMonomialBasis<dimension,StorageField> Basis;

  public:
    typedef typename Basis::DomainVector DomainVector;

    OrthonormalBasis (int maxOrder)
      : basis_(), basisEval_(0)
    {
      ONB::ONBMatrix<Topology,ComputationField> onbMatrix(maxOrder);
      coeffMatrix_.fill(onbMatrix);
      std::ofstream out("coeffs.out");
      onbMatrix.print(out);
      out << " ************ " << std::endl;
      print(out,maxOrder);
    }

    const int size (unsigned int order) const
    {
      return basis_.size(order);
    }

    template <class RangeVector>
    void evaluate ( unsigned int order,
                    const DomainVector &x,
                    std::vector< RangeVector > &values ) const
    {
      basisEval_.resize( size( order ) );
      basis_.evaluate(order,x,basisEval_);
      coeffMatrix_.mult(basisEval_,values);
    }

    template <class DomainVector,class RangeVector>
    void evaluate ( unsigned int order,
                    const DomainVector &x,
                    std::vector< RangeVector > &values ) const
    {
      DomainVector bx;
      for (int d=0; d<dimension; ++d)
        field_cast(x[d], bx[ d ]);
      evaluate(order,bx,values);
    }

    void print(std::ofstream &out,int order) {
      typedef Dune::MultiIndex<dimension> MI;
      typedef Dune::StandardMonomialBasis< dimension, MI > Basis;
      Basis basis;
      const unsigned int size = basis.sizes( order )[ order ];
      std::vector< MI > y( size );
      Dune::FieldVector< MI, dimension > x;
      for (int d=0; d<dimension; ++d)
        x[d].set(d);
      basis.evaluate( order, x, &(y[0]) );
      coeffMatrix_.print(out,y);
    }
  private:
    const Basis basis_;
    mutable std::vector<StorageField> basisEval_;
    CoeffMatrix<FieldVector<StorageField,1> > coeffMatrix_;
  };
}
#endif // DUNE_ORTHONORMALBASIS_HH
