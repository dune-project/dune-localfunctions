// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_MULTIINDEX_HH
#define DUNE_MULTIINDEX_HH

#include <dune/common/fvector.hh>

#include <dune/common/field.hh>

namespace Dune
{

  // Internal Forward Declarations
  // -----------------------------

  template< int dim >
  class MultiIndex;

  template< int dim >
  std::ostream &operator<< ( std::ostream &, const MultiIndex< dim > & );



  // MultiIndex
  // ----------

  template< int dim >
  class MultiIndex
  {
    typedef MultiIndex< dim > This;

    friend std::ostream &operator<<<> ( std::ostream &, const This & );

  public:
    static const int dimension = dim;

    MultiIndex ()
      : vecZ_( 0 ),
        vecOMZ_( 0 ),
        factor_( 1. ),
        next_( 0 )
    {}
    MultiIndex (double f)
      : vecZ_( 0 ),
        vecOMZ_( 0 ),
        factor_( f ),
        next_( 0 )
    {}

    MultiIndex ( int, const MultiIndex &other )
      : vecZ_( other.vecOMZ_ ),
        vecOMZ_( other.vecZ_ ),
        factor_( other.factor_ )
    {
      const This *o = &other;
      if (o->next_)
      {
        next_ = new MultiIndex( *(other.next_) );
      }
      else
        next_ = 0;
    }

    MultiIndex ( const This &other )
      : vecZ_( other.vecZ_ ),
        vecOMZ_( other.vecOMZ_ ),
        factor_( other.factor_ )
    {
      if (other.next_)
      {
        next_ = new MultiIndex( *(other.next_) );
      }
      else
        next_ = 0;
    }

    int z(int i) const
    {
      return vecZ_[i];
    }
    int omz(int i) const
    {
      return vecOMZ_[i];
    }
    double factor() const
    {
      return factor_;
    }

    This &operator= ( const This &other )
    {
      vecZ_   = other.vecZ_;
      vecOMZ_ = other.vecOMZ_;
      factor_ = other.factor_;
      assert(!next_);
      if (other.next_)
      {
        next_ = new MultiIndex;
        next_ = other.next_;
      }
      return *this;
    }
    This &operator= ( const double f )
    {
      factor_ = f;
      return *this;
    }

    This &operator*= ( const double f )
    {
      factor_ *= f;
      return *this;
    }
    This &operator/= ( const double f )
    {
      factor_ /= f;
      return *this;
    }

    This &operator*= ( const This &other )
    {
      vecZ_   += other.vecZ_;
      vecOMZ_ += other.vecOMZ_;
      factor_ *= other.factor_;
      return *this;
    }
    This &operator/= ( const This &other )
    {
      vecZ_   -= other.vecZ_;
      vecOMZ_ -= other.vecOMZ_;
      factor_ /= other.factor_;
      return *this;
    }

    This &operator+= ( const This &other )
    {
      if (!sameMultiIndex(other))
      {
        if (next_)
          (*next_)+=other;
        else
        {
          next_ = new This(other);
        }
      }
      else
        factor_ += other.factor_;
      return *this;
    }
    This &operator-= ( const This &other )
    {
      if (!sameMultiIndex(other))
      {
        if (next_)
          next_+=other;
        else
        {
          next_ = new This(other);
        }
      }
      else
        factor_ -= other.factor_;
      return *this;
    }

    This operator* ( const double f ) const
    {
      This z = *this;
      return (z *= f);
    }
    This operator/ ( const double f ) const
    {
      This z = *this;
      return (z /= f);
    }

    This operator* ( const This &other ) const
    {
      This z = *this;
      return (z *= other);
    }
    This operator/ ( const This &other ) const
    {
      This z = *this;
      return (z /= other);
    }

    This operator+ ( const This &other ) const
    {
      This z = *this;
      return (z += other);
    }
    This operator- ( const This &other ) const
    {
      This z = *this;
      return (z -= other);
    }

    void set ( int d, int power = 1 )
    {
      vecZ_[ d ] = power;
    }

    int absZ () const
    {
      int ret = 0;
      for( int i = 0; i < dimension; ++i )
        ret += std::abs( vecZ_[ i ] );
      return ret;
    }

    int absOMZ() const
    {
      int ret = 0;
      for( int i = 0; i < dimension; ++i )
        ret += std::abs( vecOMZ_[ i ] );
      return ret;
    }

    bool sameMultiIndex(const MultiIndex &ind)
    {
      for( int i = 0; i < dimension; ++i )
      {
        if ( vecZ_[i] != ind.vecZ_[i] ||
             vecOMZ_[i] != vecOMZ_[i] )
          return false;
      }
      return true;
    }

  private:
    typedef Dune::FieldVector< int, dimension > Vector;

    Vector vecZ_;
    Vector vecOMZ_;
    double factor_;

    This *next_;
  };

  template <int dim>
  MultiIndex<dim> operator* ( const double f,
                              const MultiIndex<dim> &m)
  {
    MultiIndex<dim> z = m;
    return (z *= f);
  }
  template <int dim>
  MultiIndex<dim> operator/ ( const double f,
                              const MultiIndex<dim> &m)
  {
    MultiIndex<dim> z = m;
    return (z /= f);
  }

  template <int d>
  std::ostream &operator<<(std::ostream& out,const std::vector<MultiIndex<d> >& y) {
    for (unsigned int r=0; r<y.size(); ++r) {
      out << "f_" << r << "(" << char('a');
      for (int i=1; i<d; ++i)
        out << "," << char('a'+i);
      out << ")=";
      out << y[r] << std::endl;
    }
    return out;
  }
  template <int d,int dimR>
  std::ostream &operator<<(std::ostream& out,
                           const std::vector<Dune::FieldVector<MultiIndex<d>,dimR> >& y) {
    for (unsigned int k=0; k<y.size(); ++k) {
      out << "f_" << k << "(" << char('a');
      for (int i=1; i<d; ++i)
        out << "," << char('a'+i);
      out << ") = ( ";
      out << y[k][0] ;
      for (unsigned int r=1; r<dimR; ++r) {
        out << " , " << y[k][r] ;
      }
      out << " )" << std::endl;
    }
    return out;
  }
  template <int d>
  std::ostream &operator<<(std::ostream& out,const MultiIndex<d>& mi) {
    if (mi.next_)
    {
      assert( &mi != mi.next_ );
      out << *(mi.next_) << " + ";
    }
    if (mi.absZ()==0 && std::abs(mi.factor())<1e-10)
      out << "0";
    else if (mi.absZ()==0)
      out << mi.factor();
    else {
      if ( std::abs(mi.factor()-1.)>1e-10)
        out << mi.factor();
      int absVal = 0;
      for (int i=0; i<d; ++i) {
        if (mi.vecZ_[i]==0)
          continue;
        else if (mi.vecZ_[i]==1)
          out << char('a'+i);
        /*
           else if (mi.vecZ_[i]>0)
           out << char('a'+i) << "**(" << mi.vecZ_[i] << ")";
           else if (mi.vecZ_[i]<0)
           out << char('a'+i) << "**(" << mi.vecZ_[i] << ")";
           absVal += mi.vecZ_[i];
           if (absVal<mi.absZ()) out << "*";
         */
        else if (mi.vecZ_[i]>0)
          out << char('a'+i) << "^" << mi.vecZ_[i] << "";
        else if (mi.vecZ_[i]<0)
          out << char('a'+i) << "^" << mi.vecZ_[i] << "";
        absVal += mi.vecZ_[i];
        if (absVal<mi.absZ()) out << "";
      }
    }
    if (mi.absOMZ()>0) {
      for (int i=0; i<=mi.absZ(); ++i) {
        if (mi.vecOMZ_[i]==0)
          continue;
        else if (mi.vecOMZ_[i]==1)
          out << (1-char('a'+i));
        else if (mi.vecOMZ_[i]>0)
          out << (1-char('a'+i)) << "^" << mi.vecOMZ_[i];
        else if (mi.vecOMZ_[i]<0)
          out << (1-char('a'+i)) << "^" << mi.vecOMZ_[i];
        if (i==mi.absZ()+1) out << "*";
      }
    }
    return out;
  }

  template< int dim >
  struct Unity< MultiIndex< dim > >
  {
    typedef MultiIndex< dim > Field;

    operator Field () const
    {
      return Field();
    }

    Field operator- ( const Field &other ) const
    {
      return Field( 1, other );
    }

    Field operator/ ( const Field &other ) const
    {
      return Field() / other;
    }
  };



  template< int dim >
  struct Zero< MultiIndex< dim > >
  {
    typedef MultiIndex< dim > Field;

    // zero does not acutally exist
    operator Field ()
    {
      assert( false );
      return Field();
    }
  };

  template< int dim >
  bool operator< ( const Zero< MultiIndex< dim > > &, const MultiIndex< dim > & )
  {
    return true;
  }

  template< int dim >
  bool operator< ( const MultiIndex< dim > &f, const Zero< MultiIndex< dim > > & )
  {
    return true;
  }

}

#endif // #ifndef DUNE_MULTIINDEX_HH
