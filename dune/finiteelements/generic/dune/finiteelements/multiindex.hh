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
      assert(!other.next_);
      if (other.next_)
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

    ~MultiIndex()
    {
      remove();
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
      remove();
      vecZ_   = other.vecZ_;
      vecOMZ_ = other.vecOMZ_;
      factor_ = other.factor_;
      if (other.next_)
        next_ = new MultiIndex(*(other.next_));
      return *this;
    }
    This &operator= ( const double f )
    {
      remove();
      vecZ_ = 0;
      vecOMZ_ = 0;
      factor_ = f;
      return *this;
    }

    This &operator*= ( const double f )
    {
      factor_ *= f;
      if (next_)
        (*next_) *= f;
      return *this;
    }
    This &operator/= ( const double f )
    {
      factor_ /= f;
      if (next_)
        (*next_) /= f;
      return *this;
    }

    This &operator*= ( const This &other )
    {
      assert(!other.next_);
      vecZ_   += other.vecZ_;
      vecOMZ_ += other.vecOMZ_;
      factor_ *= other.factor_;
      if (next_)
        (*next_) *= other;
      return *this;
    }
    This &operator/= ( const This &other )
    {
      assert(!other.next_);
      vecZ_   -= other.vecZ_;
      vecOMZ_ -= other.vecOMZ_;
      factor_ /= other.factor_;
      if (next_)
        (*next_) /= other;
      return *this;
    }

    This &operator+= ( const This &other )
    {
      assert(!other.next_);
      if (std::abs(other.factor_)<1e-10)
        return *this;
      if (std::abs(factor_)<1e-10)
      {
        *this = other;
        return *this;
      }
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
      assert(!other.next_);
      if (!sameMultiIndex(other))
      {
        if (next_)
          next_-=other;
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
    void remove()
    {
      if (next_)
      {
        next_->remove();
        delete next_;
        next_ = 0;
      }
    }

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
      out << "f_{" << r << "}(" << char('a');
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
    out << "\\begin{eqnarray*}" << std::endl;
    for (unsigned int k=0; k<y.size(); ++k) {
      out << "f_{" << k << "}(" << char('a');
      for (int i=1; i<d; ++i)
        out << "," << char('a'+i);
      out << ") &=& ( ";
      out << y[k][0] ;
      for (unsigned int r=1; r<dimR; ++r) {
        out << " , " << y[k][r] ;
      }
      out << " ) \\\\" << std::endl;
    }
    out << "\\end{eqnarray*}" << std::endl;
    return out;
  }
  template <int d,int dimR1,int dimR2>
  std::ostream &operator<<(std::ostream& out,
                           const std::vector<Dune::FieldMatrix<MultiIndex<d>,dimR1,dimR2> >& y) {
    out << "\\begin{eqnarray*}" << std::endl;
    for (unsigned int k=0; k<y.size(); ++k) {
      for (int q=0; q<dimR2; q++) {
        out << "d_{" << char('a'+q) << "}f_{" << k << "}(" << char('a');
        for (int i=1; i<d; ++i)
          out << "," << char('a'+i);
        out << ") &=& ( ";
        out << y[k][0][q] ;
        for (unsigned int r=1; r<dimR1; ++r) {
          out << " , " << y[k][r][q] ;
        }
        out << " ) \\\\" << std::endl;
      }
    }
    out << "\\end{eqnarray*}" << std::endl;
    return out;
  }
  template <int d>
  std::ostream &operator<<(std::ostream& out,const MultiIndex<d>& val)
  {
    bool first = true;
    const MultiIndex<d> *m = &val;
    do {
      if (m->absZ()==0 && std::abs(m->factor())<1e-10)
      {
        if (!m->next_ || !first)
        {
          out << "0";
          break;
        }
        else {
          m = m->next_;
          continue;
        }
      }
      if (m->factor()>0 && !first)
        out << " + ";
      else if (m->factor()<0)
        if (!first)
          out << " - ";
        else
          out << "- ";
      else
        out << "  ";
      first = false;
      double f = std::abs(m->factor());
      if (m->absZ()==0)
        out << f;
      else {
        if ( std::abs(f)<1e-10)
          out << 0;
        else
        {
          if ( std::abs(f-1.)>1e-10)
            out << f;
          int absVal = 0;
          for (int i=0; i<d; ++i) {
            if (m->vecZ_[i]==0)
              continue;
            else if (m->vecZ_[i]==1)
              out << char('a'+i);
            else if (m->vecZ_[i]>0)
              out << char('a'+i) << "^" << m->vecZ_[i] << "";
            else if (m->vecZ_[i]<0)
              out << char('a'+i) << "^" << m->vecZ_[i] << "";
            absVal += m->vecZ_[i];
            if (absVal<m->absZ()) out << "";
          }
        }
      }
      /*
         if (mi.absOMZ()>0) {
         for (int i=0;i<=mi.absZ();++i) {
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
       */
      m = m->next_;
    } while (m);
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
      return Field(0);
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
