// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ORTHONORMALBASIS_HH
#define DUNE_ORTHONORMALBASIS_HH
const unsigned int Precision = 1024;
#include "../coeffmatrix.hh"
#include "../monomialbasis.hh"
#include "orthonormalcompute.hh"
namespace Dune
{
  template <class Topology,int maxOrder,class F>
  struct ONBMatrix {
    enum {dim = Topology::dimension};
    typedef amp::ampf<Precision> scalar_t;
    typedef ap::template_1d_array< scalar_t > vec_t;
    typedef ap::template_2d_array< scalar_t > mat_t;
    ONBMatrix()
      : calc(1)
    {
      MultiIndex<dim> geo;
      calc.compute(geo);
    }
    int colSize(int row) const {
      return row;
    }
    int rowSize() const {
      return calc.res.gethighbound(1)-1;
    }
    void set(int r,int c,double v) const {
      v = calc.res(r,c).toDouble();
    }
    CalcCoeffs<dim,maxOrder> calc;
  };
  template <class Topology,int maxOrder,class F>
  class OrthonormalBasis
  {
    typedef OrthonormalBasis<Topology,maxOrder,F> This;
    typedef StandardMonomialBasis<Topology::dimension,F> Basis;

  public:
    typedef F Field;

    typedef typename Basis::DomainVector DomainVector;
    typedef typename Basis::RangeVector RangeVector;

    OrthonormalBasis ()
      : basis_(), basisEval_(0)
    {
      ONBMatrix<Topology,maxOrder,Field> onbMatrix;
      coeffMatrix_.fill(onbMatrix);
    }

    const int size (unsigned int order) const
    {
      return basis_.size(order);
    }

    void evaluate ( unsigned int order,
                    const DomainVector &x,
                    std::vector< RangeVector > &values ) const
    {
      assert(order<=maxOrder);
      basisEval_.resize(size(order));
      basis_.evaluate(order,x,basisEval_);
      mult(basisEval_,values);
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
