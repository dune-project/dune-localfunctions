// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_MONOMIALBASIS_HH
#define DUNE_MONOMIALBASIS_HH

#include <dune/common/fvector.hh>

#include <dune/grid/genericgeometry/topologytypes.hh>

namespace Dune
{

  // MonomialBasis
  // -------------

  template< class Topology, class F >
  class MonomialBasis;

  template< class F >
  class MonomialBasis< GenericGeometry::Point, F >
  {
    typedef MonomialBasis< GenericGeometry::Point, F > This;

  public:
    typedef GenericGeometry::Point Topology;
    typedef F Field;

    static const unsigned int dimDomain = Topology::dimension;
    static const unsigned int dimRange = 1;

    typedef FieldVector< Field, dimDomain > DomainVector;
    typedef FieldVector< Field, dimRange > RangeVector;

  private:
    mutable int maxOrder_;
    mutable int *sizes_;     // sizes_[k] is number of basis functions of exactly order k
                             // sizes_[0] = 1,....
    mutable int *numBaseFunctions_; // nBF[k] = sizes_[0]+...+sizes_[k] is number of
                                    // basis functions up to order k

    template< int dimD >
    void evaluate ( const unsigned int order,
                    const FieldVector< Field, dimD > &x,
                    const unsigned int *const offsets,
                    RangeVector *const values ) const
    {
      values[0] = Field(1);
    }

    int maxOrder() const {
      return maxOrder_;
    }
    void computeSizes(int order) const {
      maxOrder_ = order;
      delete [] sizes_;
      delete [] numBaseFunctions_;
      int *sizes_            = new int [order+1];
      int *numBaseFunctions_ = new int [order+1];
      newSizes[0] = 1;
      numBaseFunctions_[0] = 1;
      for (int k=1; k<=order; ++k) {
        sizes_[k]              = 0;
        numBaseFunctions_[k]   = 1;
      }
    }
  public:
    const int* sizes(int order) const {
      if ( order>maxOrder() ) {
        computeSizes(order);
      }
      return numBaseFunctions_;
    }
    void evaluate ( const unsigned int order,
                    const DomainVector &x,
                    RangeVector *const values ) const
    {
      evaluate( order, x, sizes(order), values );
    }
  };

  template< class BaseTopology, class F >
  class MonomialBasis< GenericGeometry::Prism< BaseTopology >, F >
  {
    typedef MonomialBasis< GenericGeometry::Prism< BaseTopology >, F > This;

  public:
    typedef GenericGeometry::Prism< BaseTopology > Topology;
    typedef F Field;

    static const unsigned int dimDomain = Topology::dimension;
    static const unsigned int dimRange = 1;

    typedef FieldVector< Field, dimDomain > DomainVector;
    typedef FieldVector< Field, dimRange > RangeVector;

  private:
    MonomialBasis< BaseTopology, Field > baseBasis_;
    int *sizes_;     // sizes_[k] is number of basis functions of exactly order k
                     // sizes_[0] = 1,....
    int *numBaseFunctions_; // nBF[k] = sizes_[0]+...+sizes_[k] is number of
                            // basis functions up to order k

  public:
    MonomialBasis () : sizes_(0), numBaseFunctions_(0)
    {
      computeSizes(2);
    }

  private:
    template< int dimD >
    void evaluate ( const unsigned int order,
                    const FieldVector< Field, dimD > &x,
                    const unsigned int *const offsets,
                    RangeVector *const values ) const
    {
      const Field& z = x[dimDomain-1];

      baseBasis_.evaluate( order, x, offsets, values ); // fill first column
      const int baseSizes = baseBasis.sizes_;

      RangeVector *row0 = values;
      for( int k = 1; k <= order; ++k )
      {
        RangeVector * row1 = values+offsets[k-1];
        RangeVector * it = row1+basesSizes[k];
        RangeVector * colkEnd = row1+(k+1)*basesSizes[k];
        for( ; it!=colkEnd; ++row1,++it ) {
          *it = z * (*row1);
        }
        RangeVector *const row1End = row1+sizes_[k];
        for( ; it!=row1End; ++row0,++it ) {
          *it = z * (*row0);
        }
        row0 = row1;
      }
    }
    int maxOrder() const {
      return baseBasis_.maxOrder();
    }
    void computeSizes(int order) const {
      delete [] sizes_;
      delete [] numBaseFunctions_;
      int *sizes_            = new int [order+1];
      int *numBaseFunctions_ = new int [order+1];
      baseBasis_.computeSizes(order);
      const int baseSizes = baseBasis.sizes_;
      const int baseNBF   = baseBasis.numBaseFunctions_;
      newSizes[0] = 1;
      numBaseFunctions_[0] = 1;
      for (int k=1; k<=order; ++k) {
        sizes_[k]              = baseNBF[k]+k*basesSizes[k];
        numBaseFunctions_[k]   = numBaseFunctions_[k-1]+baseNBF[k];
      }
    }
  public:
    const int* sizes(int order) const {
      if ( order>maxOrder() ) {
        computeSizes(order);
      }
      return numBaseFunctions_;
    }
    void evaluate ( const unsigned int order,
                    const DomainVector &x,
                    RangeVector *const values ) const
    {
      evaluate( order, x, sizes(), values );
    }
  };

  template< class BaseTopology, class F >
  class MonomialBasis< GenericGeometry::Pyramid< BaseTopology >, F >
  {
    typedef MonomialBasis< GenericGeometry::Pyramid< BaseTopology >, F > This;

  public:
    typedef GenericGeometry::Pyramid< BaseTopology > Topology;
    typedef F Field;

    static const unsigned int dimDomain = Topology::dimension;
    static const unsigned int dimRange = 1;

    typedef FieldVector< Field, dimDomain > DomainVector;
    typedef FieldVector< Field, dimRange > RangeVector;

  private:
    MonomialBasis< BaseTopology, Field > baseBasis_;
    mutable int *sizes_;     // sizes_[k] is number of basis functions of exactly order k
                             // sizes_[0] = 1,....
    mutable int *numBaseFunctions_; // nBF[k] = sizes_[0]+...+sizes_[k] is number of
                                    // basis functions up to order k
  public:
    MonomialBasis () : sizes_(0), numBaseFunctions_(0)
    {
      computeSizes(2);
    }

  private:
    template< int dimD >
    void evaluateSimplex ( const int order,
                           const FieldVector< Field, dimD > &x,
                           const int *const offsets,
                           RangeVector *const values ) const
    {
      const Field& z = x[dimDomain-1];

      baseBasis_.evaluate( order, x, offsets, values ); // fill first column
      const int baseSizes = baseBasis.sizes_;

      RangeVector *row0 = values;
      for( int k = 1; k <= order; ++k )
      {
        RangeVector *const row1 = values+offsets[k-1];
        RangeVector *const row1End = row1+sizes_[k];
        for( RangeVector *it = row1+baseSizes[k]; it!=row1End; ++row0,++it ) {
          *it = z * (*row0);
        }
        row0 = row1;
      }
    }

    template< int dimD >
    void evaluatePyramid ( const unsigned int order,
                           const FieldVector< Field, dimD > &x,
                           const unsigned int *const offsets,
                           RangeVector *const values ) const
    {
      const Field& z = x[dimDomain-1];
      const Field& omz = Field(1)-z;
      const Field& invomz = Field(1)/omz;
      const FieldVector< Field, dimDomain-1 > y;
      for (int i=0; i<dimDomain-1; i++)
        y[i] = x[i] * invomz;

      baseBasis_.evaluate( order, y, offsets, values ); // fill first column
      const int baseSizes = baseBasis.sizes_;

      RangeVector *row0 = values;
      const Field& omzk = omz;
      for( int k = 1; k <= order; ++k )
      {
        RangeVector *const row1 = values+offsets[k-1];
        RangeVector *const row1End = row1+sizes_[k];
        RangeVector *const col0End = row1+baseSizes[k];
        RangeVector *it = row1;
        for( ; it!=col0End; ++it ) {
          *it = *it * omzk;
        }
        for( ; it!=row1End; ++row0,++it ) {
          *it = z * (*row0);
        }
        row0 = row1;
        omzk *= omz;
      }
    }

    template< int dimD >
    void evaluate ( const unsigned int order,
                    const FieldVector< Field, dimD > &x,
                    const unsigned int *const offsets,
                    RangeVector *const values ) const
    {
      if( GenericGeometry::IsSimplex< Topology >::value )
        evaluateSimplex( order, x, offsets, values );
      else
        evaluatePyramid( order, x, offsets, values );
    }
    int maxOrder() const {
      return baseBasis_.maxOrder();
    }
    void computeSizes(int order) const {
      delete [] sizes_;
      delete [] numBaseFunctions_;
      int *sizes_            = new int [order+1];
      int *numBaseFunctions_ = new int [order+1];
      baseBasis_.computeSizes(order);
      const int baseNBF = baseBasis.numBaseFunctions_;
      newSizes[0] = 1;
      numBaseFunctions_[0] = 1;
      for (int k=1; k<=order; ++k) {
        sizes_[k]              = baseNBF[k];
        numBaseFunctions_[k]   = numBaseFunctions_[k-1]+baseNBF[k];
      }
    }
  public:
    const int* sizes(int order) const {
      if ( order>maxOrder() ) {
        computeSizes(order);
      }
      return numBaseFunctions_;
    }
    void evaluate ( const unsigned int order,
                    const DomainVector &x,
                    RangeVector *const values ) const
    {
      evaluate( order, x, sizes(order), values );
    }
  };
}

#endif
