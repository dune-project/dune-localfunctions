// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <iostream>

#include <dune/finiteelements/lagrangepoints.hh>
#include <dune/finiteelements/monomialbasis.hh>

#ifndef TOPOLOGY
#error "TOPOLOGY not defined."
#endif

using namespace Dune::GenericGeometry;

namespace Dune {
  template <int dim>
  struct MultiIndex {
    MultiIndex( ) {
      vecZ_   = Vector(0);
      vecOMZ_ = Vector(0);
    }
    void set(int d) {
      vecZ_[d] = 1;
    }
    MultiIndex(int ,const MultiIndex& other)
      : vecZ_(other.vecOMZ_), vecOMZ_(other.vecZ_)
    {}
    MultiIndex(const MultiIndex& other)
      : vecZ_(other.vecZ_), vecOMZ_(other.vecOMZ_)
    {}
    MultiIndex& operator=(const MultiIndex& other) {
      vecZ_   = other.vecZ_;
      vecOMZ_ = other.vecOMZ_;
      return *this;
    }
    MultiIndex& operator*=(const MultiIndex& other) {
      vecZ_   += other.vecZ_;
      vecOMZ_ += other.vecOMZ_;
      return *this;
    }
    MultiIndex& operator/=(const MultiIndex& other) {
      vecZ_   -= other.vecZ_;
      vecOMZ_ -= other.vecOMZ_;
      return *this;
    }
    MultiIndex operator* (const MultiIndex& b) const {
      MultiIndex z = *this;
      return (z*=b);
    }
    MultiIndex operator/ (const MultiIndex& b) const {
      MultiIndex z = *this;
      return (z/=b);
    }
    int absZ() const {
      int ret=0;
      for (int i=0; i<dim; ++i)
        ret+=vecZ_[i];
      return ret;
    }
    int absOMZ() const {
      int ret=0;
      for (int i=0; i<dim; ++i)
        ret+=std::abs(vecOMZ_[i]);
      return ret;
    }
    typedef Dune::FieldVector<int,dim> Vector;
    Vector vecZ_;
    Vector vecOMZ_;
  };
  template <int d>
  std::ostream &operator<<(std::ostream& out,const MultiIndex<d>& mi) {
    for (int i=0; i<d; ++i) {
      if (mi.absZ()==0)
        out << "1   ";
      else if (mi.vecZ_[i]==0)
        out << "    ";
      else if (mi.vecZ_[i]==1)
        out << char('a'+i) << "   ";
      else if (mi.vecZ_[i]>0)
        out << char('a'+i) << "^" << mi.vecZ_[i] << " ";
      else if (mi.vecZ_[i]<0)
        out << char('a'+i) << "^" << mi.vecZ_[i] << " ";
    }
    out << "   ";
    for (int i=0; i<d; ++i)
      out << "(1-" << char('a'+i) << ")^" << mi.vecOMZ_[i] << " ";
    return out;
  }

  template <int d>
  struct Unity<MultiIndex<d> > {
    typedef MultiIndex<d> Field;
    operator Field () const
    {
      return Field();
    }
    Field operator- (const Field& other) const {
      return Field(1,other);
    }
    Field operator/ (const Field& other) const {
      return Field()/other;
    }
  };
  template <int d>
  struct Zero<MultiIndex<d> > {
    typedef MultiIndex<d> Field;
    operator Field() { // does not acctualy exist
      assert(0);
      return Field();
    }
  };
  template <int d>
  bool operator<(const Zero< MultiIndex<d> >& ,const MultiIndex<d>& f) {
    return true;
  }
  template <int d>
  bool operator<(const MultiIndex<d>& f, const Zero< MultiIndex<d> >&) {
    return true;
  }
}

int main ( int argc, char **argv )
{
  if( argc < 2 )
  {
    std::cerr << "Usage: " << argv[ 0 ] << " <p>" << std::endl;
    return 1;
  }

  int p = atoi( argv[ 1 ] );

  typedef TOPOLOGY Topology;

  Dune::MonomialBasis< Topology, double > basis;
  const unsigned int size = basis.sizes( p )[ p ];
  std::vector< Dune::FieldVector< double, 1 > > y( size );

  typedef Dune::LagrangePoints< Topology, double > LagrangePoints;
  LagrangePoints points( p );

  std::cout << "Number of base functions:  " << size << std::endl;
  std::cout << "Number of Lagrange points: " << points.size() << std::endl;

  const LagrangePoints::iterator end = points.end();
  for( LagrangePoints::iterator it = points.begin(); it != end; ++it )
  {
    basis.evaluate( p, *it, y );
    std::cout << "x = " << *it << ":" << std::endl;
    for( unsigned int i = 0; i < size; ++i )
      std::cout << "    y[ " << i << " ] = " << y[ i ] << std::endl;
  }
  {
    typedef Dune::StandardBiMonomialBasis< 3,double > Basis;
    Basis basis;
    const unsigned int size = basis.sizes( p )[ p ];
    std::vector< Dune::FieldVector< double, 1 > > y( size );

    typedef Dune::LagrangePoints< Basis::Topology, double > LagrangePoints;
    LagrangePoints points( p );

    const LagrangePoints::iterator end = points.end();
    for( LagrangePoints::iterator it = points.begin(); it != end; ++it )
    {
      basis.evaluate( p, *it, y );
      std::cout << "x = " << *it << ":" << std::endl;
      for( unsigned int i = 0; i < size; ++i )
        std::cout << "    y[ " << i << " ] = " << y[ i ] << std::endl;
    }
  }
  {
    enum {dim=Topology::dimension};
    typedef Dune::MultiIndex<dim> MI;
    typedef Dune::MonomialBasis< Topology, MI  > Basis;
    Basis basis;
    const unsigned int size = basis.sizes( p )[ p ];
    std::vector< Dune::FieldVector< MI,1> > y( size );
    Dune::FieldVector< MI, dim > x;
    for (int d=0; d<dim; ++d)
      x[d].set(d);
    basis.evaluate( p, x, y );
    for (int i=0; i<y.size(); ++i)
      std::cout << y[i] << std::endl;
  }
}
