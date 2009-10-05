// NOTE: The current revision of this file was left untouched when the DUNE source files were reindented!
// NOTE: It contained invalid syntax that could not be processed by uncrustify.

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
    template< int dimD >
    void evaluate ( const unsigned int order,
                    const FieldVector< Field, dimD > &x,
                    const unsigned int *const offsets,
                    double *const values ) const
    {}
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
    MonomialBasis< BaseTopology, Field > baseBasis;
    int *sizes_;

  private:
    template< int dimD >
    void evaluate ( const unsigned int order,
                    const FieldVector< Field, dimD > &x,
                    const unsigned int *const offsets,
                    double *const values ) const
    {
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
    MonomialBasis< BaseTopology, Field > baseBasis;
    int *sizes_;

  public:
    MonomialBasis ()
    {}

  private:
    template< int dimD >
    void evaluateSimplex ( const int order,
                           const FieldVector< Field, dimD > &x,
                           const int *const offsets,
                           double *const values ) const
    {
      baseBasis.evaluate( order, x, offsets, values );
      const int baseSizes = baseBasis.sizes_;

      for( int k = 0; k < p; ++k )
      {
        // evaluate all polynomials of order k
        for( int j = k; j >= 0; --j )
          

      }
    }

    template< int dimD >
    void evaluatePyramid ( const unsigned int order,
                           const FieldVector< Field, dimD > &x,
                           const unsigned int *const offsets,
                           double *const values ) const
    {
    }

    template< int dimD >
    void evaluate ( const unsigned int order,
                    const FieldVector< Field, dimD > &x,
                    const unsigned int *const offsets,
                    double *const values ) const
    {
      if( GenericGeometry::IsSimplex< Topology >::value )
        evaluateSimplex( order, x, offsets, values );
      else
        evaluatePyramid( order, x, offsets, values );
    }
  };
}

#endif
