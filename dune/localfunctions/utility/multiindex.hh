// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_MULTIINDEX_HH
#define DUNE_MULTIINDEX_HH

#include <vector>
#include <ostream>

#include <dune/common/fvector.hh>

#include <dune/localfunctions/utility/field.hh>

namespace Dune
{
  /****************************************************************
  * Provide a Field class which can be used in evaluation methods
  * to produce MultiIndex presentation of polynomials.
  ****************************************************************/
  // Internal Forward Declarations
  // -----------------------------

  template< int dim, class Field >
  class MultiIndex;

  template< int dim, class Field >
  std::ostream &operator<< ( std::ostream &, const MultiIndex< dim,Field > & );



  // MultiIndex
  // ----------

  template< int dim,class Field >
  class MultiIndex
  {
    typedef MultiIndex< dim, Field > This;

    friend std::ostream &operator<<<> ( std::ostream &, const This & );

  public:
    static const int dimension = dim;

    MultiIndex ()
      : vecZ_( 0 ),
        vecOMZ_( 0 ),
        factor_( 1. ),
        next_( 0 )
    {}
    template <class F>
    explicit MultiIndex (const F &f)
      : vecZ_( 0 ),
        vecOMZ_( 0 ),
        factor_( field_cast<Field>(f) ),
        next_( 0 )
    {}

    MultiIndex ( int, const This &other )
      : vecZ_( other.vecOMZ_ ),
        vecOMZ_( other.vecZ_ ),
        factor_( other.factor_ )
    {
      assert(!other.next_);
      if (other.next_)
      {
        next_ = new This( *(other.next_) );
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
        next_ = new This( *(other.next_) );
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
    const Field &factor() const
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
        next_ = new This(*(other.next_));
      return *this;
    }
    This &operator= ( const Zero<This> &f )
    {
      remove();
      vecZ_ = 0;
      vecOMZ_ = 0;
      factor_ = 0.;
      return *this;
    }
    This &operator= ( const Unity<This> &f )
    {
      remove();
      vecZ_ = 0;
      vecOMZ_ = 0;
      factor_ = 1.;
      return *this;
    }
    template <class F>
    This &operator= ( const F &f )
    {
      remove();
      vecZ_ = 0;
      vecOMZ_ = 0;
      factor_ = field_cast<Field>(f);
      return *this;
    }

    bool operator== (const This &other) const
    {
      assert(!next_ && !other.next_);
      return (vecZ_==other.vecZ_ && vecOMZ_==other.vecOMZ_ && factor_==other.factor_);
    }

    template <class F>
    This &operator*= ( const F &f )
    {
      factor_ *= field_cast<Field>(f);
      if (next_)
        (*next_) *= f;
      return *this;
    }
    template <class F>
    This &operator/= ( const F &f )
    {
      factor_ /= field_cast<Field>(f);
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

    template <class F>
    This operator* ( const F &f ) const
    {
      This z = *this;
      return (z *= f);
    }
    template <class F>
    This operator/ ( const F &f ) const
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

    bool sameMultiIndex(const This &ind)
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
    Field factor_;

    This *next_;
  };

  template <int dim, class Field, class F>
  MultiIndex<dim,Field> operator* ( const F &f,
                                    const MultiIndex<dim,Field> &m)
  {
    MultiIndex<dim,Field> z = m;
    return (z *= f);
  }
  template <int dim, class Field, class F>
  MultiIndex<dim,Field> operator/ ( const F &f,
                                    const MultiIndex<dim,Field> &m)
  {
    MultiIndex<dim,Field> z = m;
    return (z /= f);
  }

  template <int d, class F>
  std::ostream &operator<<(std::ostream& out,const std::vector<MultiIndex<d,F> >& y) {
    for (unsigned int r=0; r<y.size(); ++r) {
      out << "f_{" << r << "}(" << char('a');
      for (int i=1; i<d; ++i)
        out << "," << char('a'+i);
      out << ")=";
      out << y[r] << std::endl;
    }
    return out;
  }
  template <int d,class F,int dimR>
  std::ostream &operator<<(std::ostream& out,
                           const std::vector<Dune::FieldVector<MultiIndex<d,F>,dimR> >& y) {
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
  template <int d,class F,int dimR1,int dimR2>
  std::ostream &operator<<(std::ostream& out,
                           const std::vector<Dune::FieldMatrix<MultiIndex<d,F>,dimR1,dimR2> >& y) {
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
  template <int d, class F>
  std::ostream &operator<<(std::ostream& out,const MultiIndex<d,F>& val)
  {
    bool first = true;
    const MultiIndex<d,F> *m = &val;
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
      F f = std::abs(m->factor());
      if (m->absZ()==0)
        out << f;
      else {
        if ( std::abs(f)<1e-10)
          out << 0;
        else
        {
          F f_1(f);
          f_1 -= 1.; // better Unity<F>();
          if ( std::abs(f_1)>1e-10)
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

  template< int dim, class F>
  struct Unity< MultiIndex< dim, F > >
  {
    typedef MultiIndex< dim, F > Field;

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



  template< int dim, class F >
  struct Zero< MultiIndex< dim,F > >
  {
    typedef MultiIndex< dim,F > Field;

    // zero does not acutally exist
    operator Field ()
    {
      return Field(0);
    }
  };

  template< int dim, class Field >
  bool operator< ( const Zero< MultiIndex< dim,Field > > &, const MultiIndex< dim,Field > & )
  {
    return true;
  }

  template< int dim, class Field >
  bool operator< ( const MultiIndex< dim, Field > &f, const Zero< MultiIndex< dim,Field > > & )
  {
    return true;
  }

}

#endif // #ifndef DUNE_MULTIINDEX_HH
