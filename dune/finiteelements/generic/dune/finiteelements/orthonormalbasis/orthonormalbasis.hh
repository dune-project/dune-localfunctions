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
  template <class Topology,class F>
  class OrthonormalBasis
  {
    enum {dim = Topology::dimension};
    typedef OrthonormalBasis<Topology,F> This;
    typedef StandardMonomialBasis<dim,F> Basis;
    static const unsigned int Precision = 1024;

  public:
    typedef F Field;

    typedef typename Basis::DomainVector DomainVector;
    typedef typename Basis::RangeVector RangeVector;

    OrthonormalBasis (int maxOrder)
      : basis_(), basisEval_(0)
    {
      ONB::ONBMatrix<Topology,amp::ampf<Precision> > onbMatrix(maxOrder);
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

    void evaluate ( unsigned int order,
                    const DomainVector &x,
                    std::vector< RangeVector > &values ) const
    {
      basisEval_.resize(size(order));
      basis_.evaluate(order,x,basisEval_);
      mult(basisEval_,values);
    }

    void print(std::ofstream &out,int order) {
      typedef Dune::MultiIndex<dim> MI;
      typedef Dune::MonomialBasis< Topology, MI  > Basis;
      Basis basis;
      const unsigned int size = basis.sizes( order )[ order ];
      std::vector< Dune::FieldVector< MI,1> > y( size );
      Dune::FieldVector< MI, dim > x;
      for (int d=0; d<dim; ++d)
        x[d].set(d);
      basis.evaluate( order, x, y );
      coeffMatrix_.print(out,y);
    }
  private:
    void mult(const std::vector< RangeVector > &x,
              std::vector< RangeVector > &y) const {
      coeffMatrix_.mult(x,y);
    }
    const Basis basis_;
    mutable std::vector<RangeVector> basisEval_;
    CoeffMatrix<Field> coeffMatrix_;
  };
}
#endif // DUNE_ORTHONORMALBASIS_HH
