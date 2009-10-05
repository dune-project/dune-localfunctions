// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_MONOMIALBASIS_HH
#define DUNE_MONOMIALBASIS_HH

#include <vector>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <dune/grid/genericgeometry/topologytypes.hh>
#include <dune/grid/genericgeometry/misc.hh>

#include <dune/common/field.hh>

#include <dune/finiteelements/basisprovider.hh>
#include <dune/finiteelements/multiindex.hh>

namespace Dune
{

  // Internal Forward Declarations
  // -----------------------------

  template< class Topology, class F >
  class MonomialBasis;



  // MonomialBasisImpl
  // -----------------

  template< class Topology, class F >
  class MonomialBasisImpl;

  template< class F >
  class MonomialBasisImpl< GenericGeometry::Point, F >
  {
    typedef MonomialBasisImpl< GenericGeometry::Point, F > This;

  public:
    typedef GenericGeometry::Point Topology;
    typedef F Field;

    static const unsigned int dimDomain = Topology::dimension;

    typedef FieldVector< Field, dimDomain > DomainVector;

  private:
    friend class MonomialBasis< Topology, Field >;
    friend class MonomialBasisImpl< GenericGeometry::Prism< Topology >, Field >;
    friend class MonomialBasisImpl< GenericGeometry::Pyramid< Topology >, Field >;

    mutable unsigned int maxOrder_;
    // sizes_[ k ]: number of basis functions of exactly order k
    mutable unsigned int *sizes_;
    // numBaseFunctions_[ k ] = sizes_[ 0 ] + ... + sizes_[ k ]
    mutable unsigned int *numBaseFunctions_;

    MonomialBasisImpl ()
      : sizes_( 0 ),
        numBaseFunctions_( 0 )
    {
      computeSizes( 2 );
    }

    ~MonomialBasisImpl ()
    {
      delete[] sizes_;
      delete[] numBaseFunctions_;
    }

    template< int dimD >
    void evaluate ( const unsigned int order,
                    const FieldVector< Field, dimD > &x,
                    const unsigned int *const offsets,
                    Field *const values ) const
    {
      values[ 0 ] = Unity< Field >();
    }

    void integral ( const unsigned int order,
                    const unsigned int *const offsets,
                    Field *const values ) const
    {
      values[ 0 ] = Unity< Field >();
    }

    unsigned int maxOrder () const
    {
      return maxOrder_;
    }

    void computeSizes ( unsigned int order ) const
    {
      maxOrder_ = order;

      delete [] sizes_;
      delete [] numBaseFunctions_;
      sizes_            = new unsigned int [ order+1 ];
      numBaseFunctions_ = new unsigned int [ order+1 ];

      sizes_[ 0 ] = 1;
      numBaseFunctions_[ 0 ] = 1;
      for( unsigned int k = 1; k <= order; ++k )
      {
        sizes_[ k ]            = 0;
        numBaseFunctions_[ k ] = 1;
      }
    }
  };

  template< class BaseTopology, class F >
  class MonomialBasisImpl< GenericGeometry::Prism< BaseTopology >, F >
  {
    typedef MonomialBasisImpl< GenericGeometry::Prism< BaseTopology >, F > This;

  public:
    typedef GenericGeometry::Prism< BaseTopology > Topology;
    typedef F Field;

    static const unsigned int dimDomain = Topology::dimension;

    typedef FieldVector< Field, dimDomain > DomainVector;

  private:
    friend class MonomialBasis< Topology, Field >;
    friend class MonomialBasisImpl< GenericGeometry::Prism< Topology >, Field >;
    friend class MonomialBasisImpl< GenericGeometry::Pyramid< Topology >, Field >;

    MonomialBasisImpl< BaseTopology, Field > baseBasis_;
    // sizes_[ k ]: number of basis functions of exactly order k
    mutable unsigned int *sizes_;
    // numBaseFunctions_[ k ] = sizes_[ 0 ] + ... + sizes_[ k ]
    mutable unsigned int *numBaseFunctions_;

    MonomialBasisImpl ()
      : sizes_( 0 ),
        numBaseFunctions_( 0 )
    {
      computeSizes( 2 );
    }

    ~MonomialBasisImpl ()
    {
      delete[] sizes_;
      delete[] numBaseFunctions_;
    }

    template< int dimD >
    void evaluate ( const unsigned int order,
                    const FieldVector< Field, dimD > &x,
                    const unsigned int *const offsets,
                    Field *const values ) const
    {
      const Field &z = x[ dimDomain-1 ];

      // fill first column
      baseBasis_.evaluate( order, x, offsets, values );
      const unsigned int *const baseSizes = baseBasis_.sizes_;

      Field *row0 = values;
      for( unsigned int k = 1; k <= order; ++k )
      {
        Field *const row1begin = values + offsets[ k-1 ];
        Field *const colkEnd = row1begin + (k+1)*baseSizes[ k ];
        assert( (unsigned int)(colkEnd - values) <= offsets[ k ] );
        Field *const row1End = row1begin + sizes_[ k ];
        assert( (unsigned int)(row1End - values) <= offsets[ k ] );

        Field *row1 = row1begin;
        Field *it;
        for( it = row1begin + baseSizes[ k ]; it != colkEnd; ++row1, ++it )
          *it = z * (*row1);
        for( ; it != row1End; ++row0, ++it )
          *it = z * (*row0);
        row0 = row1;
      }
    }

    void integral ( const unsigned int order,
                    const unsigned int *const offsets,
                    Field *const values ) const
    {
      /*
         // fill first column
         baseBasis_.integral( order, offsets, values );
         const unsigned int *const baseSizes = baseBasis_.sizes_;

         Field *row0 = values;
         for( unsigned int k = 1; k <= order; ++k )
         {
         Field *const row1begin = values + offsets[ k-1 ];
         Field *const row1End = row1begin + sizes_[ k ];
         assert( (unsigned int)(row1End - values) <= offsets[ k ] );

         Field *row1 = row1begin;
         Field *it = row1begin + baseSizes[ k ];
         for( unsigned int j = 1; j <= k; ++j )
         {
          Field *const end = it + baseSizes[ k ];
          assert( (unsigned int)(end - values) <= offsets[ k ] );
          for( ; it != end; ++row1, ++it )
         *it = (Field( j ) / Field( j+1 )) * (*row1);
         }
         for( ; it != row1End; ++row0, ++it )
         *it = (Field( k ) / Field( k+1 )) * (*row0);
         row0 = row1;
         }
       */
    }

    unsigned int maxOrder() const
    {
      return baseBasis_.maxOrder();
    }

    void computeSizes ( unsigned int order ) const
    {
      delete[] sizes_;
      delete[] numBaseFunctions_;
      sizes_            = new unsigned int[ order+1 ];
      numBaseFunctions_ = new unsigned int[ order+1 ];

      baseBasis_.computeSizes( order );
      const unsigned int *const baseSizes = baseBasis_.sizes_;
      const unsigned int *const baseNBF   = baseBasis_.numBaseFunctions_;

      sizes_[ 0 ] = 1;
      numBaseFunctions_[ 0 ] = 1;
      for( unsigned int k = 1; k <= order; ++k )
      {
        sizes_[ k ]            = baseNBF[ k ] + k*baseSizes[ k ];
        numBaseFunctions_[ k ] = numBaseFunctions_[ k-1 ] + sizes_[ k ];
      }
    }
  };

  template< class BaseTopology, class F >
  class MonomialBasisImpl< GenericGeometry::Pyramid< BaseTopology >, F >
  {
    typedef MonomialBasisImpl< GenericGeometry::Pyramid< BaseTopology >, F > This;

  public:
    typedef GenericGeometry::Pyramid< BaseTopology > Topology;
    typedef F Field;

    static const unsigned int dimDomain = Topology::dimension;

    typedef FieldVector< Field, dimDomain > DomainVector;

  private:
    friend class MonomialBasis< Topology, Field >;
    friend class MonomialBasisImpl< GenericGeometry::Prism< Topology >, Field >;
    friend class MonomialBasisImpl< GenericGeometry::Pyramid< Topology >, Field >;

    MonomialBasisImpl< BaseTopology, Field > baseBasis_;
    // sizes_[ k ]: number of basis functions of exactly order k
    mutable unsigned int *sizes_;
    // numBaseFunctions_[ k ] = sizes_[ 0 ] + ... + sizes_[ k ]
    mutable unsigned int *numBaseFunctions_;

    MonomialBasisImpl ()
      : sizes_( 0 ),
        numBaseFunctions_( 0 )
    {
      computeSizes( 2 );
    }

    ~MonomialBasisImpl ()
    {
      delete[] sizes_;
      delete[] numBaseFunctions_;
    }

    template< int dimD >
    void evaluateSimplex ( const unsigned int order,
                           const FieldVector< Field, dimD > &x,
                           const unsigned int *const offsets,
                           Field *const values ) const
    {
      const Field &z = x[ dimDomain-1 ];

      // fill first column
      baseBasis_.evaluate( order, x, offsets, values );

      const unsigned int *const baseSizes = baseBasis_.sizes_;
      Field *row0 = values;
      for( unsigned int k = 1; k <= order; ++k )
      {
        Field *const row1 = values+offsets[ k-1 ];
        Field *const row1End = row1+sizes_[ k ];
        assert( (unsigned int)(row1End - values) <= offsets[ k ] );
        for( Field *it = row1 + baseSizes[ k ]; it != row1End; ++row0, ++it )
          *it = z * (*row0);
        row0 = row1;
      }
    }

    template< int dimD >
    void evaluatePyramid ( const unsigned int order,
                           const FieldVector< Field, dimD > &x,
                           const unsigned int *const offsets,
                           Field *const values ) const
    {
      const Field &z = x[ dimDomain-1 ];
      Field omz = Unity< Field >() - z;

      if( Zero<Field>() < omz )
      {
        const Field invomz = Unity<Field>() / omz;
        FieldVector< Field, dimDomain-1 > y;
        for( unsigned int i = 0; i < dimDomain-1; ++i )
          y[ i ] = x[ i ] * invomz;

        // fill first column
        baseBasis_.evaluate( order, y, offsets, values );

        const unsigned int *const baseSizes = baseBasis_.sizes_;
        Field *row0 = values;
        Field omzk = omz;
        for( unsigned int k = 1; k <= order; ++k )
        {
          Field *const row1 = values + offsets[ k-1 ];
          Field *const row1End = row1 + sizes_[ k ];
          assert( (unsigned int)(row1End - values) <= offsets[ k ] );
          Field *const col0End = row1 + baseSizes[ k ];
          Field *it = row1;
          for( ; it != col0End; ++it )
            *it = (*it) * omzk;
          for( ; it != row1End; ++row0, ++it )
            *it = z * (*row0);
          row0 = row1;
          omzk *= omz;
        }
      }
      else {
        Field *it = values;
        for( unsigned int k = 0; k <= order; ++k )
        {
          Field *const rowEnd = it + (sizes_[ k ] - 1);
          for( ; it != rowEnd; ++it )
            *it = Zero<Field>();
          *it = Unity<Field>();
          ++it;
        }
      }

    }

    template< int dimD >
    void evaluate ( const unsigned int order,
                    const FieldVector< Field, dimD > &x,
                    const unsigned int *const offsets,
                    Field *const values ) const
    {
      if( GenericGeometry::IsSimplex< Topology >::value )
        evaluateSimplex( order, x, offsets, values );
      else
        evaluatePyramid( order, x, offsets, values );
    }

    void integral ( const unsigned int order,
                    const unsigned int *const offsets,
                    Field *const values ) const
    {
      /*
         // fill first column
         baseBasis_.integral( order, offsets, values );

         const unsigned int *const baseSizes = baseBasis_.sizes_;

         Field *const col0End = values + baseSizes[ 0 ];
         for( Field *it = values; it != col0End; ++it )
         *it *= Field( 1 ) /  Field( double(dimDomain) );  // ??? double cast due to error in Linker
         Field *row0 = values;

         for( unsigned int k = 1; k <= order; ++k )
         {
          const Field factor = (Field( 1 ) / Field( k + dimDomain ));

          Field *const row1 = values+offsets[ k-1 ];
          Field *const col0End = row1 + baseSizes[ k ];
          Field *it = row1;
          for( ; it != col0End; ++it )
         *it *= factor;
          for( unsigned int i = 1; i <= k; ++i )
          {
            Field *const end = it + baseSizes[ k-i ];
            assert( (unsigned int)(end - values) <= offsets[ k ] );
            for( ; it != end; ++row0, ++it )
         *it = (factor * Field( i )) * (*row0);
          }
          row0 = row1;
         }
       */
    }

    unsigned int maxOrder() const
    {
      return baseBasis_.maxOrder();
    }

    void computeSizes ( unsigned int order ) const
    {
      delete[] sizes_;
      delete[] numBaseFunctions_;
      sizes_            = new unsigned int[ order+1 ];
      numBaseFunctions_ = new unsigned int[ order+1 ];

      baseBasis_.computeSizes( order );

      const unsigned int *const baseNBF = baseBasis_.numBaseFunctions_;
      sizes_[ 0 ] = 1;
      numBaseFunctions_[ 0 ] = 1;
      for( unsigned int k = 1; k <= order; ++k )
      {
        sizes_[ k ]            = baseNBF[ k ];
        numBaseFunctions_[ k ] = numBaseFunctions_[ k-1 ] + sizes_[ k ];
      }
    }
  };



  // MonomialBasis
  // -------------

  template< class Topology, class F >
  class MonomialBasis
    : public MonomialBasisImpl< Topology, F >
  {
    typedef MonomialBasis< Topology, F > This;
    typedef MonomialBasisImpl< Topology, F > Base;

  public:
    static const unsigned int dimension = Base::dimDomain;

    typedef typename Base::Field Field;

    typedef typename Base::DomainVector DomainVector;

    MonomialBasis ()
      : Base()
    {}

    const unsigned int *sizes ( unsigned int order ) const
    {
      if( order > Base::maxOrder() )
        Base::computeSizes( order );
      return Base::numBaseFunctions_;
    }

    const unsigned int size ( unsigned int order ) const
    {
      return sizes( order )[ order ];
    }

    void evaluate ( const unsigned int order,
                    const DomainVector &x,
                    Field *const values ) const
    {
      Base::evaluate( order, x, sizes( order ), values );
    }

    void evaluate ( const unsigned int order,
                    const DomainVector &x,
                    FieldVector<Field,1> *const values ) const
    {
      evaluate( order, x, reinterpret_cast< Field * >( values ) );
    }

    template <class RangeVector>
    void evaluate ( const unsigned int order,
                    const DomainVector &x,
                    std::vector< RangeVector > &values ) const
    {
      evaluate( order, x, &(values[ 0 ]) );
    }

    void integral ( const unsigned int order,
                    Field *const values ) const
    {
      Base::integral( order, sizes( order ), values );
    }

    void integral ( const unsigned int order,
                    FieldVector<Field,1> *const values ) const
    {
      integral( order, reinterpret_cast< Field * >( values ) );
    }

    template <class RangeVector>
    void integral ( const unsigned int order,
                    std::vector< RangeVector > &values ) const
    {
      integral( order, &(values[ 0 ]) );
    }
  private:
    MonomialBasis(const This&);
    This& operator=(const This&);
  };
  template< int dim, class F >
  class VirtualMonomialBasis
  {
    typedef VirtualMonomialBasis< dim, F > This;

  public:
    typedef F Field;
    static const int dimension = dim;

    typedef FieldVector<Field,dimension> DomainVector;

    virtual ~VirtualMonomialBasis() {}

    virtual const unsigned int *sizes ( unsigned int order ) const = 0;

    const unsigned int size ( unsigned int order ) const
    {
      return sizes( order )[ order ];
    }

    virtual void evaluate ( const unsigned int order,
                            const DomainVector &x,
                            Field *const values ) const = 0;
    void evaluate ( const unsigned int order,
                    const DomainVector &x,
                    FieldVector<Field,1> *const values ) const
    {
      evaluate( order, x, reinterpret_cast< Field * >( values ) );
    }
    template <class RangeVector>
    void evaluate ( const unsigned int order,
                    const DomainVector &x,
                    std::vector< RangeVector > &values ) const
    {
      evaluate( order, x, &(values[ 0 ]) );
    }

    virtual void integral ( const unsigned int order,
                            Field *const values ) const = 0;
    void integral ( const unsigned int order,
                    FieldVector<Field,1> *const values ) const
    {
      integral( order, reinterpret_cast< Field * >( values ) );
    }
    template <class RangeVector>
    void integral ( const unsigned int order,
                    std::vector< RangeVector > &values ) const
    {
      integral( order, &(values[ 0 ]) );
    }
  };
  template< class Topology, class F >
  class VirtualMonomialBasisImpl
    : public VirtualMonomialBasis< Topology::dimension, F >
  {
    typedef VirtualMonomialBasis< Topology::dimension, F > Base;
    typedef VirtualMonomialBasisImpl< Topology, F > This;

  public:
    typedef typename Base::Field Field;
    typedef typename Base::DomainVector DomainVector;

    const unsigned int *sizes ( unsigned int order ) const
    {
      return basis_.sizes(order);
    }

    void evaluate ( const unsigned int order,
                    const DomainVector &x,
                    Field *const values ) const
    {
      basis_.evaluate(order,x,values);
    }

    void integral ( const unsigned int order,
                    Field *const values ) const
    {
      basis_.integral(order,values);
    }

  private:
    MonomialBasis<Topology,Field> basis_;
  };

  template< int dim, class F >
  struct MonomialBasisCreator
  {
    typedef F StorageField;
    typedef VirtualMonomialBasis<dim,StorageField> Basis;
    static const int dimension = dim;
    typedef unsigned int Key;

    template <class Topology>
    struct Maker
    {
      static void apply(unsigned int order,Basis* &basis)
      {
        basis = new VirtualMonomialBasisImpl<Topology,StorageField>;
      }
    };
  };

  template< int dim, class SF >
  struct MonomialBasisProvider
    : public BasisProvider<MonomialBasisCreator<dim,SF> >
  {};

  // StdMonomialTopology
  // -------------------

  template <int dim>
  struct StdMonomialTopology {
    typedef StdMonomialTopology<dim-1> BaseType;
    typedef GenericGeometry::Pyramid<typename BaseType::Type> Type;
  };
  template <>
  struct StdMonomialTopology<0> {
    typedef GenericGeometry::Point Type;
  };
  template <int dim>
  struct StdBiMonomialTopology {
    typedef GenericGeometry::Prism<typename StdBiMonomialTopology<dim-1>::Type> Type;
  };
  template <>
  struct StdBiMonomialTopology<0> {
    typedef GenericGeometry::Point Type;
  };
  template< int dim,class F >
  class StandardMonomialBasis
    : public MonomialBasis< typename StdMonomialTopology<dim>::Type, F >
  {
  public:
    typedef typename StdMonomialTopology<dim>::Type Topology;
    static const int dimension = dim;
  private:
    typedef StandardMonomialBasis< dim, F > This;
    typedef MonomialBasis< Topology, F > Base;
  public:
    StandardMonomialBasis ()
      : Base()
    {}
  };
  template< int dim, class F >
  class StandardBiMonomialBasis
    : public MonomialBasis< typename StdBiMonomialTopology<dim>::Type, F >
  {
  public:
    typedef typename StdBiMonomialTopology<dim>::Type Topology;
    static const int dimension = dim;
  private:
    typedef StandardBiMonomialBasis< dim, F > This;
    typedef MonomialBasis< Topology, F > Base;
  public:
    StandardBiMonomialBasis ()
      : Base()
    {}
  };

  template <class F,int dimD,int dimR,unsigned int deriv>
  struct Tensor
  {
    typedef Tensor<F,dimD,dimR,deriv-1> Base;
    static const int single = Base::single*dimD;
    static const int all = Base::all+single;
    typedef Dune::FieldVector<Dune::FieldVector<F,single>,dimR> Single;
    typedef Dune::FieldVector<Dune::FieldVector<F,all>,dimR> All;
  };
  template <class F,int dimD,int dimR>
  struct Tensor<F,dimD,dimR,2>
  {
    struct HJV
    {
      typedef FieldVector<FieldMatrix<F,dimD,dimD>,dimR> Hessian;
      typedef FieldMatrix<F,dimR,dimD> Jacobian;
      typedef FieldVector<F,dimR> Value;
      Hessian hessian;
      Jacobian jacobian;
      Value value;
    };
    typedef Tensor<F,dimD,dimR,1> Base;
    static const int single = Base::single*dimD;
    static const int all = Base::all+single;
    typedef typename HJV::Hessian Single;
    typedef HJV All;
  };
  template <class F,int dimD,int dimR>
  struct Tensor<F,dimD,dimR,1>
  {
    struct JV
    {
      typedef FieldMatrix<F,dimR,dimD> Jacobian;
      typedef FieldVector<F,dimR> Value;
      Jacobian jacobian;
      Value value;
    };
    typedef Tensor<F,dimD,dimR,0> Base;
    static const int single = Base::single*dimD;
    static const int all = Base::all+single;
    typedef typename JV::Jacobian Single;
    typedef JV All;
  };
  template <class F,int dimD,int dimR>
  struct Tensor<F,dimD,dimR,0>
  {
    static const int single = 1;
    static const int all = 1;
    typedef Dune::FieldVector<F,dimR> Single;
    typedef Dune::FieldVector<F,dimR> All;
  };

  template <class B>
  struct StandardEvaluator
  {
    typedef B Basis;
    typedef typename Basis::Field Field;
    typedef typename Basis::DomainVector DomainVector;
    typedef std::vector<Field> Container;
    static const int dimension = Basis::dimension;
    template <class BlockType>
    struct BaseIterator
    {
      static const int blockSize = sizeof(BlockType)/sizeof(Field);
      typedef BlockType RangeVector;
      typedef typename Container::const_iterator CIter;
      BaseIterator(const Container &container)
        : pos_(container.begin()), end_(container.end())
      {}
      const RangeVector &operator*() const
      {
        assert(!done());
        return reinterpret_cast<const RangeVector&>(*pos_);
      }
      const RangeVector *operator->() const
      {
        assert(!done());
        return reinterpret_cast<const RangeVector*>(pos_);
      }
      bool done() const
      {
        return pos_==end_;
      }
      BaseIterator &operator++()
      {
        assert(blockSize == 1);
        pos_ += blockSize;
        return *this;
      }
      BaseIterator &operator+=(unsigned int skip)
      {
        assert(blockSize == 1);
        pos_ += skip*blockSize;
        return *this;
      }
    private:
      CIter pos_;
      const CIter end_;
    };

    StandardEvaluator(const Basis &basis,unsigned int order)
      : basis_(basis),
        order_(order),
        container_(basis.size(order))
    {}
    template <unsigned int deriv>
    struct Iterator
    {
      typedef BaseIterator<typename Tensor<Field,dimension,1,deriv>::All> All;
      typedef BaseIterator<typename Tensor<Field,dimension,1,deriv>::Single> Single;
    };
    typename Iterator<0>::Single evaluate(const DomainVector &x)
    {
      basis_.evaluate(order_,x,container_);
      return typename Iterator<0>::Single(container_);
    }
    template <unsigned int deriv>
    typename Iterator<deriv>::Single evaluate(const DomainVector &x)
    {
      basis_.evaluate(order_,x,container_);
      return typename Iterator<deriv>::Single(container_);
    }
    typename Iterator<1>::Single jacobian(const DomainVector &x)
    {
      basis_.evaluate(order_,x,container_);
      return typename Iterator<1>::Single(container_);
    }
    template <unsigned int deriv>
    typename Iterator<deriv>::All evaluateAll(const DomainVector &x)
    {
      basis_.evaluateAll<deriv>(order_,x,container_);
      return typename Iterator<deriv>::All(container_);
    }
    unsigned int order() const
    {
      return order_;
    }
    unsigned int size() const
    {
      return basis_.size(order);
    }
  private:
    StandardEvaluator(const StandardEvaluator&);
    const Basis &basis_;
    unsigned int order_;
    Container container_;
  };

  template <class B,class F>
  struct MultiIndexEvaluator
  {
    typedef B Basis;
    typedef F Field;
    typedef std::vector<typename Basis::Field> Container;
    static const int dimension = Basis::dimension;
    typedef Dune::FieldVector<Field,dimension> DomainVector;
    template <class BlockType,int deriv>
    struct BaseIterator
    {
      static const int blockSize = sizeof(BlockType)/sizeof(Field);
      typedef BlockType RangeVector;
      typedef typename Container::const_iterator CIter;
      BaseIterator(const std::vector<DomainVector> &x,
                   const Container &container)
        : pos_(container.begin()), end_(container.end()),
          x_(x)
      {
        set();
      }
      const RangeVector &operator*() const
      {
        assert(!done());
        return reinterpret_cast<const RangeVector&>(val_);
      }
      const RangeVector * const operator->() const
      {
        assert(!done());
        return reinterpret_cast<const RangeVector* const>(&val_);
      }
      bool done() const
      {
        return pos_==end_;
      }
      BaseIterator operator++()
      {
        pos_ += 1;
        if (!done())
          set();
        return *this;
      }
      BaseIterator &operator+=(unsigned int skip)
      {
        pos_ += skip;
        if (!done())
          set();
        return *this;
      }
    private:
      void set()
      {
        if (deriv==0)
        {
          val_[0]=1;
          for (int d=0; d<dimension; ++d)
          {
            unsigned int o = pos_->z(d);
            assert( o<x_.size() );
            val_[0]  *= x_[o][d];
          }
        }
        else if (deriv==1)
        {
          for (int i=0; i<dimension; ++i)
          {
            unsigned int o = pos_->z(i);
            if ( o == 0)
              val_[i] = 0.;
            else
            {
              val_[i] = o;
              for (int d=0; d<dimension; ++d)
              {
                unsigned int o = pos_->z(d);
                o -= (d==i);
                assert( o<x_.size() );
                val_[i]  *= x_[o][d];
              }
            }
          }
          if (blockSize>dimension || deriv==0)
          {
            val_[dimension]=1;
            for (int d=0; d<dimension; ++d)
            {
              unsigned int o = pos_->z(d);
              assert( o<x_.size() );
              val_[dimension]  *= x_[o][d];
            }
          }
        }
      }
      CIter pos_;
      const CIter end_;
      const std::vector<DomainVector> &x_;
      Field val_[blockSize];
    };

    MultiIndexEvaluator(const Basis &basis,unsigned int order)
      : basis_(basis),
        order_(order),
        container_(basis.size(order)),
        x_(order+1)
    {
      typename Basis::DomainVector x;
      for( int i = 0; i < dimension; ++i )
        x[ i ].set( i, 1 );
      basis_.evaluate(order_,x,container_);
    }
    template <unsigned int deriv>
    struct Iterator
    {
      typedef BaseIterator<typename Tensor<Field,dimension,1,deriv>::All,deriv> All;
      typedef BaseIterator<typename Tensor<Field,dimension,1,deriv>::Single,deriv> Single;
    };
    typename Iterator<0>::Single evaluate(const DomainVector &x)
    {
      setX(x);
      return typename Iterator<0>::Single(x_,container_);
    }
    template <unsigned int deriv>
    typename Iterator<deriv>::Single evaluate(const DomainVector &x)
    {
      setX(x);
      return typename Iterator<deriv>::Single(x_,container_);
    }
    template <unsigned int deriv>
    typename Iterator<deriv>::All evaluateAll(const DomainVector &x)
    {
      setX(x);
      return typename Iterator<deriv>::All(x_,container_);
    }
    typename Iterator<1>::Single jacobian(const DomainVector &x)
    {
      setX(x);
      return typename Iterator<1>::Single(x_,container_);
    }
    void setX(const DomainVector &x)
    {
      for (int d=0; d<dimension; ++d)
      {
        x_[0][d] = 1;
        for (unsigned int i=1; i<=order_; ++i) {
          x_[i][d]=x_[i-1][d]*x[d];
        }
      }
    }
    unsigned int order() const
    {
      return order_;
    }
    unsigned int size() const
    {
      return basis_.size(order);
    }
  private:
    MultiIndexEvaluator(const MultiIndexEvaluator&);
    const Basis &basis_;
    unsigned int order_;
    Container container_;
    std::vector<DomainVector> x_;
  };
}

#endif
