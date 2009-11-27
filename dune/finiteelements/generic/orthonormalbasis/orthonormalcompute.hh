// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ORTHONORMALCOMPUTE_HH
#define DUNE_ORTHONORMALCOMPUTE_HH

#include <cassert>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <map>

#include <dune/common/fmatrix.hh>

#include <dune/finiteelements/generic/math/matrix.hh>

#include <dune/grid/genericgeometry/topologytypes.hh>

#include <dune/finiteelements/generic/common/monomialbasis.hh>
#include <dune/finiteelements/generic/common/multiindex.hh>

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
    static int compute(const Dune::MultiIndex<dim, scalar_t>& alpha,
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
    static int compute(const Dune::MultiIndex<dim,scalar_t>& alpha,
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
    static int compute(const Dune::MultiIndex<dim,scalar_t>& alpha,
                       scalar_t& p,scalar_t& q) {
      p = 1.;
      q = 1.;
      return 0;
    }
  };
  template <class Topology, class scalar_t>
  struct ONBMatrix : Dune::LFEMatrix<scalar_t>
  {
    typedef Dune::LFEMatrix<scalar_t> Base;

    typedef std::vector< scalar_t > vec_t;
    typedef Dune::LFEMatrix< scalar_t > mat_t;
    static const unsigned int dim=Topology::dimension;

    explicit ONBMatrix( const unsigned int order )
    {
      // get all multiindecies for monomial basis
      typedef Dune::MultiIndex<dim,scalar_t> MI;
      typedef Dune::StandardMonomialBasis< dim, MI  > Basis;
      Basis basis( order );
      const unsigned int size = basis.size( );
      std::vector< Dune::FieldVector< MI,1> > y( size );
      Dune::FieldVector< MI, dim > x;
      for (unsigned int i=0; i<dim; ++i) {
        x[i].set(i);
      }
      basis.evaluate( x, y );
      // set bounds of data
      Base::resize(size,size);
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
      gramSchmidt();
    }
    template <class Vector>
    void row( const unsigned int row, Vector &vec ) const
    {
      // transposed matrix is required
      assert(row<Base::cols());
      for (unsigned int i=0; i<Base::rows(); ++i)
        field_cast(Base::operator()(i,row), vec[i]);
    }
    bool test() {
      bool ret = true;
      const unsigned int N = Base::rows();
      // Nun schauen wir noch, ob wirklich C^T S C = E gilt
      scalar_t prod;
      for (unsigned int i=0; i<N; ++i) {
        for (unsigned int j=0; j<N; j++) {
          prod = 0;
          for (unsigned int k=0; k<N; k++) {
            for (unsigned int l=0; l<N; l++) {
              prod += Base::operator()(l,i)*S(l,k)*Base::operator()(k,j);
            }
          }
          assert((i==j) ? 1 : fabs( field_cast<double>(prod) )<1e-10);
        }
      }
      return ret;
    }
  private:
    void sprod(int col1,int col2, scalar_t& ret)
    {
      ret = 0;
      for (int k=0; k<=col1; k++) {
        for (int l=0; l<=col2; l++) {
          ret += Base::operator()(l,col2)*this->S(l,k)*Base::operator()(k,col1);
        }
      }
    }
    void vmul(unsigned int col, unsigned int rowEnd, scalar_t &ret)
    {
      for (unsigned int i=0; i<=rowEnd; ++i)
        Base::operator()(i,col) *= ret;
    }
    void vsub(unsigned int coldest, unsigned int colsrc,
              unsigned int rowEnd, scalar_t &ret)
    {
      for (unsigned int i=0; i<=rowEnd; ++i)
        Base::operator()(i,coldest) -= ret*Base::operator()(i,colsrc);
    }
    void gramSchmidt()
    {
      const unsigned int N = Base::rows();
      scalar_t s;
      for (unsigned int i=0; i<N; ++i)
        for (unsigned int j=0; j<N; ++j)
          Base::operator()(i,j) = (i==j) ? 1 : 0;
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

    vec_t d;
    mat_t S;
  };

}
#endif
