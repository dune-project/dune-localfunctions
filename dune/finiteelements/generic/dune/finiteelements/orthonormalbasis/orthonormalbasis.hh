// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ORTHONORMALBASIS_HH
#define DUNE_ORTHONORMALBASIS_HH
#include "../coeffmatrix.hh"
#include "../monomialbasis.hh"
const unsigned int Precision = 1024;
#include "orthonormalcompute.hh"
namespace Dune
{
  template <class Topology>
  struct TopologyToMultiIndex;
  template <class Base>
  struct TopologyToMultiIndex<GenericGeometry::Pyramid<Base> > {
    typedef TopologyToMultiIndex<Base> BaseType;
    template <int dim>
    static void set(MultiIndex<dim>& mi) {
      BaseType::set(mi);
      mi.set(Base::dimension,1);
    }
  };
  template <class Base>
  struct TopologyToMultiIndex<GenericGeometry::Prism<Base> > {
    typedef TopologyToMultiIndex<Base> BaseType;
    template <int dim>
    static void set(MultiIndex<dim>& mi) {
      BaseType::set(mi);
      mi.set(Base::dimension,2);
    }
  };
  template <>
  struct TopologyToMultiIndex<GenericGeometry::Pyramid<GenericGeometry::Point> > {
    template <int dim>
    static void set(MultiIndex<dim>& mi) {
      mi.set(0,1);
    }
  };
  template <>
  struct TopologyToMultiIndex<GenericGeometry::Prism<GenericGeometry::Point> > {
    template <int dim>
    static void set(MultiIndex<dim>& mi) {
      mi.set(0,2);
    }
  };
  // **********************************************
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
      TopologyToMultiIndex<Topology>::set(geo);
      calc.compute(geo);
    }
    int colSize(int row) const {
      return row;
    }
    int rowSize() const {
      return calc.res.gethighbound(1)-1;
    }
    void set(int r,int c,double v) const {
      v = calc.res(r+1,c+1).toDouble();
    }
    void set(int r,int c,std::string v) const {
      v = amp::ampf<128>(calc.res(r+1,c+1)).toDec();
    }
    void print(std::ostream& out) {
      int N = rowSize();
      for (int i=0; i<=N; ++i) {
        out << "Polynomial : " << i << std::endl;
        for (int j=0; j<=colSize(i); j++) {
          double v = 0;
          set(i,j,v);
          if (fabs(v)<1e-20)
            out << 0 << "\t\t" << std::flush;
          else {
            std::string v;
            set(i,j,v);
            out << v << "\t\t" << std::flush;
          }
        }
        for (int j=i+1; j<=N; j++) {
          assert(fabs(res(j,i).toDouble())<1e-10);
        }
        out << std::endl;
      }
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
      std::ofstream out("coeffs.out");
      onbMatrix.print(out);
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
