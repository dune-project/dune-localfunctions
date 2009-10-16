// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ORTHONORMALCOMPUTE_HH
#define DUNE_ORTHONORMALCOMPUTE_HH

#include <cassert>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <map>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <dune/alglib/vector.hh>
#include <dune/alglib/matrix.hh>

#include <dune/grid/genericgeometry/topologytypes.hh>

#include <dune/finiteelements/generic/monomialbasis.hh>
#include <dune/finiteelements/multiindex.hh>

namespace ONBCompute
{

  template <class scalar_t>
  scalar_t factorial(int start,int end) {
    scalar_t ret(1);
    for (int j=start; j<=end; j++)
      ret*=j;
    return ret;
  }
  template <class Topology>
  struct Integral;
  template <class Base>
  struct Integral<Dune::GenericGeometry::Pyramid<Base> > {
    enum {d = Base::dimension+1};
    template <int dim,class scalar_t>
    static int compute(const Dune::MultiIndex<dim>& alpha,
                       scalar_t& p,scalar_t& q) {
      int i = alpha.z(d-1);
      int ord = Integral<Base>::compute(alpha,p,q);
      p *= factorial<scalar_t>(1,i);
      q *= factorial<scalar_t>(d+ord,d+ord+i);
      return ord+i;
    }
  };
  template <class Base>
  struct Integral<Dune::GenericGeometry::Prism<Base> > {
    enum {d = Base::dimension+1};
    template <int dim,class scalar_t>
    static int compute(const Dune::MultiIndex<dim>& alpha,
                       scalar_t& p,scalar_t& q) {
      int i = alpha.z(d-1);
      int ord = Integral<Base>::compute(alpha,p,q);
      Integral<Base>::compute(alpha,p,q);
      p *= 1.;
      q *= (i+1.);
      return ord+i;
    }
  };
  template <>
  struct Integral<Dune::GenericGeometry::Point> {
    template <int dim,class scalar_t>
    static int compute(const Dune::MultiIndex<dim>& alpha,
                       scalar_t& p,scalar_t& q) {
      p = 1.;
      q = 1.;
      return 0;
    }
  };
  template <class scalar_t>
  struct Compute
  {
    typedef Dune::AlgLib::Vector< scalar_t > vec_t;
    typedef Dune::AlgLib::Matrix< scalar_t > mat_t;


    /**************************************************/
    template <class Topology>
    struct CalcCoeffsBase {
      enum {dim=Topology::dimension};
      vec_t d;
      mat_t S;
      mat_t& res;
      CalcCoeffsBase(mat_t& pres) : res(pres) {}
      virtual ~CalcCoeffsBase() {}
      bool test() {
        bool ret = true;
        const unsigned int N = this->res.rows();
        // Nun schauen wir noch, ob wirklich C^T S C = E gilt
        scalar_t prod;
        for (unsigned int i=0; i<N; ++i) {
          for (unsigned int j=0; j<N; j++) {
            prod = 0;
            for (unsigned int k=0; k<N; k++) {
              for (unsigned int l=0; l<N; l++) {
                prod += res(l,i)*S(l,k)*res(k,j);
              }
            }
            assert((i==j) ? 1 : fabs( field_cast<double>(prod) )<1e-10);
          }
        }
        return ret;
      }
      virtual void compute() = 0;
      void compute(int ord) {
        // get all multiindecies for monomial basis
        typedef Dune::MultiIndex<dim> MI;
        typedef Dune::StandardMonomialBasis< dim, MI  > Basis;
        Basis basis( ord );
        const unsigned int size = basis.size( );
        std::vector< Dune::FieldVector< MI,1> > y( size );
        Dune::FieldVector< MI, dim > x;
        for (int i=0; i<dim; ++i) {
          x[i].set(i);
        }
        basis.evaluate( x, y );
        // set bounds of data
        res.resize(size,size);
        S.resize(size,size);
        d.resize(size);
        // Aufstellen der Matrix fuer die Bilinearform xSy: S_ij = int_A x^(i+j)
        scalar_t p,q;
        for( unsigned int i=0; i<size; ++i )
        {
          for( unsigned int j=0; j<size; j++ )
          {
            Integral<Topology>::compute(y[i]*y[j],p,q);
            S(i,j) = p;
            S(i,j) /= q;
          }
        }
        compute();
        assert(test());
        return;
      }
    };
    template <class Topology>
    struct CalcCoeffsGM : public CalcCoeffsBase<Topology> {
      CalcCoeffsGM(mat_t& pres)
        : CalcCoeffsBase<Topology>(pres) {}
      virtual ~CalcCoeffsGM() {}
      void sprod(int col1,int col2, scalar_t& ret) {
        ret = 0;
        for (int k=0; k<=col1; k++) {
          for (int l=0; l<=col2; l++) {
            ret += this->res(l,col2)*this->S(l,k)*this->res(k,col1);
          }
        }
      }
      void vmul(unsigned int col, unsigned int rowEnd, scalar_t &ret)
      {
        for (unsigned int i=0; i<=rowEnd; ++i)
          this->res(i,col) *= ret;
      }
      void vsub(unsigned int coldest, unsigned int colsrc,
                unsigned int rowEnd, scalar_t &ret)
      {
        for (unsigned int i=0; i<=rowEnd; ++i)
          this->res(i,coldest) -= ret*this->res(i,colsrc);
      }
      virtual void compute() {
        const unsigned int N = this->res.rows();
        scalar_t s;
        for (unsigned int i=0; i<N; ++i)
          for (unsigned int j=0; j<N; ++j)
            this->res(i,j) = (i==j) ? 1 : 0;
        sprod(0,0,s);
        s = 1./sqrt(s);
        vmul(0,0,s);
        for (unsigned int i=0; i<N-1; i++) {
          for (unsigned int k=0; k<=i; ++k) {
            sprod(i+1,k,s);
            vsub(i+1,k,i+1,s);
          }
          sprod(i+1,i+1,s);
          s = 1./sqrt(s);
          vmul(i+1,i+1,s);
        }
      }
    };
  public:
    template <class Topology>
    struct CalcCoeffs {
      mat_t res;
      CalcCoeffsBase<Topology> *method;
      CalcCoeffs(int use_method) {
        if (use_method==1)
          method = new CalcCoeffsGM<Topology>(res);
        else
          method = 0; // new CalcCoeffsQR<dim,ord>(res);
      }
      ~CalcCoeffs() {
        delete method;
      }
      void compute(int order) {
        method->compute(order);
      }
    };
  };

  template< class Topology, class scalar_t >
  struct ONBMatrix
  {
    static const unsigned int dimension = Topology::dimension;

    typedef ONBCompute::Compute< scalar_t > Compute;
    typedef typename Compute::vec_t vec_t;
    typedef typename Compute::mat_t mat_t;

    explicit ONBMatrix( const unsigned int maxOrder )
      : calc( 1 )
    {
      calc.compute( maxOrder );
    }

    unsigned int colSize ( const unsigned int row ) const
    {
      return row+1;
    }

    unsigned int rowSize () const
    {
      return calc.res.rows(); // calc.res.gethighbound( 1 )+1;
    }

    scalar_t operator() ( const unsigned int r, const unsigned int c ) const
    {
      return calc.res( c, r );
    }

    void print( std::ostream &out, const unsigned int N = rowSize() ) const
    {
      for( unsigned int i = 0; i < N; ++i )
      {
        out << "Polynomial : " << i << std::endl;
        for( unsigned int j = 0; j <colSize( i ); ++j )
        {
          double v = calc.res(j,i).toDouble();
          if (fabs(v)<1e-20)
            out << 0 << "\t\t" << std::flush;
          else {
            scalar_t v = calc.res(j,i);
            out << v << "\t\t" << std::flush;
          }
        }
        for( unsigned int j = colSize( i ); j < N; ++j )
          assert(fabs(calc.res(j,i).toDouble())<1e-10);
        out << std::endl;
      }
    }
    typename Compute::template CalcCoeffs<Topology> calc;
  };
}
#endif
