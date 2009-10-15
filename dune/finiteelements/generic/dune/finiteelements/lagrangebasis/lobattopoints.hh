// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LOBATTOBASIS_HH
#define DUNE_LOBATTOBASIS_HH

#include <fstream>
#include <dune/alglib/multiprecision.hh>
#include <dune/alglib/matrix.hh>
#include <dune/common/field.hh>
#include <dune/common/forloop.hh>

#include <dune/finiteelements/lagrangebasis/lagrangepoints.hh>
#include <dune/finiteelements/lagrangebasis/interpolation.hh>
#include <dune/finiteelements/generic/basisprovider.hh>
#include <dune/finiteelements/basisprint.hh>
#include <dune/finiteelements/generic/polynomialbasis.hh>

#include <dune/grid/genericgeometry/referenceelements.hh>
#include <dune/grid/genericgeometry/referencemappings.hh>
#include <dune/finiteelements/quadrature/lobattoquadrature.hh>
#include <dune/finiteelements/lagrangebasis/lagrangebasis.hh>

namespace Dune
{

  template< class Field >
  struct LobattoPoints
  {
    LobattoPoints ( unsigned int order )
      : points_( 0 )
    {
      if( order < 2 )
        return;
      points_.resize(order-1);
      typedef amp::ampf< Precision< Field >::value > MPField;
      GenericGeometry::LobattoPoints<MPField> lobatto(order+1);
      for (unsigned int i=1; i<order; ++i) {
        points_[i-1] = field_cast<Field>(lobatto.point(i));
      }
      /*
         for (unsigned int i=1;i<order;++i)
         points_[i-1] = Field(i)/Field(order);
       */
    }

    const unsigned int size()
    {
      return points_.size();
    }
    const unsigned int order()
    {
      return points_.size()+1;
    }
    const Field &point(int i)
    {
      return points_[i];
    }
    std::vector<Field> points_;
  };

  template <class Field,class Topology>
  struct LobattoInnerPoints;

  template <class Field>
  struct LobattoInnerPoints<Field,GenericGeometry::Point>
  {
    static const unsigned int dimension = 0;
    static unsigned int size(const unsigned int order)
    {
      return 1;
    }
    template <unsigned int dim>
    static unsigned int setupSimplex(
      const unsigned int iup,
      const unsigned int dimStart,
      unsigned int startPos,
      const std::vector<Field> &points1D,
      LagrangePoint< Field, dim > *points )
    {
      const unsigned int order = points1D.size()+1;
      unsigned int i = order+1-iup;
      if (i-2>=points1D.size())
      {
        return startPos;
      }
      for ( unsigned int d=0; d<dimStart; ++d )
      {
        points[startPos].point_[d] = -points1D[i-2];
      }

      return startPos+1;
    }
    template <unsigned int dim>
    static void setup(const std::vector<Field> &points1D,
                      LagrangePoint< Field, dim > *points )
    {
      points->point_[0] = Zero<Field>();
    }
  };
  template <class Field,class Base>
  struct LobattoInnerPoints<Field,GenericGeometry::Pyramid<Base> >
  {
    typedef LobattoInnerPoints<Field,Base> LobattoBase;
    static const unsigned int dimension = Base::dimension+1;
    static unsigned int size(const unsigned int order)
    {
      if (order<=dimension) return 0;
      unsigned int size=0;
      for (unsigned int o=0; o<order-dimension; ++o)
        size += LobattoBase::size(o+dimension);
      return size;
    }
    template <unsigned int dim>
    static unsigned int setupSimplex(
      const unsigned int iup,
      const unsigned int dimStart,
      unsigned int startPos,
      const std::vector<Field> &points1D,
      LagrangePoint< Field, dim > *points )
    {
      const unsigned int order = points1D.size()+1;
      unsigned int endPos = startPos;
      for (unsigned int i=2; i<=order-iup; ++i)
      {
        endPos = LobattoBase::template setupSimplex<dim>(iup+i-1,dimStart,startPos,points1D,points);
        for (unsigned int j=startPos; j<endPos; ++j)
        {
          for (unsigned int d=0; d<dimStart; ++d)
          {
            if ( d==dimStart-dimension )
              points[j].point_[d] += dimStart*points1D[i-2];
            else
              points[j].point_[d] -= points1D[i-2];
          }
        }
        startPos = endPos;
      }
      return endPos;
    }
    template <unsigned int dim>
    static void setup(const std::vector<Field> &points1D,
                      LagrangePoint< Field, dim > *points )
    {
      const unsigned int order = points1D.size()+1;
      unsigned int startPos=0,endPos=0;
      for (unsigned int i=2; i<=order; ++i)
      {
        endPos = LobattoBase::template setupSimplex<dim>(i-1,dimension,startPos,points1D,points);

        for (unsigned int j=startPos; j<endPos; ++j)
        {
          for (unsigned int d=0; d<dimension; ++d)
          {
            if ( d==0 )
              points[j].point_[d] += 1.+int(dimension)*points1D[i-2];
            else
              points[j].point_[d] += 1.-points1D[i-2];
            points[j].point_[d] /= (dimension+1);
          }
        }
        startPos = endPos;
      }
    }
  };
  template <class Field,class Base>
  struct LobattoInnerPoints<Field,GenericGeometry::Prism<Base> >
  {
    typedef LobattoInnerPoints<Field,Base> LobattoBase;
    static const unsigned int dimension = Base::dimension+1;
    static unsigned int size(const unsigned int order)
    {
      return LobattoBase::size(order)*(order-1);
    }
    template <unsigned int dim>
    static unsigned int setupSimplex(
      const unsigned int iup,
      const unsigned int dimStart,
      unsigned int startPos,
      const std::vector<Field> &points1D,
      LagrangePoint< Field, dim > *points )
    {
      const unsigned int order = points1D.size()+1;
      unsigned int endPos = startPos;
      for (unsigned int i=2; i<=order-iup; ++i)
      {
        endPos = LobattoBase::template setupSimplex<dim>(iup+i-1,dimStart,startPos,points1D,points);
        for (unsigned int j=startPos; j<endPos; ++j)
        {
          for (unsigned int d=0; d<dimStart; ++d)
          {
            if ( d==dimStart-dimension )
              points[j].point_[d] += dimStart*points1D[i-2];
            else
              points[j].point_[d] -= points1D[i-2];
          }
        }
        startPos = endPos;
      }
      return endPos;
    }
    template <unsigned int dim>
    static void setup(const std::vector<Field> &points1D,
                      LagrangePoint< Field, dim > *points )
    {
      const unsigned int order = points1D.size()+1;
      assert(dim>=dimension);
      assert(points1D.size()==order-1);
      LobattoBase::template setup<dim>(points1D,points);
      const unsigned int baseSize = LobattoBase::size(order);
      for (unsigned int q=0; q<points1D.size(); q++)
      {
        for (unsigned int i=0; i<baseSize; ++i)
        {
          const unsigned int pos = q*baseSize+i;
          for (unsigned int d=0; d<dimension-1; ++d)
            points[pos].point_[d] = points[i].point_[d];
          points[pos].point_[dimension-1]=points1D[q];
        }
      }
    }
  };

  template <class Field,unsigned int dim>
  struct LobattoPointsCreator
  {
  public:
    static const unsigned int dimension = dim;

    typedef Dune::LagrangePoints< Field, dimension > LagrangePoints;
    typedef LagrangePoints LocalCoefficients;

    typedef unsigned int Key;

    template <class Topology>
    struct Setup
    {
      template <int pdim>
      struct InitCodim
      {
        static const unsigned int codim = dimension-pdim;
        static void apply(const unsigned int order,
                          const std::vector<Field> &points1D,
                          LagrangePoints &points)
        {
          const unsigned int size = GenericGeometry::Size<Topology,codim>::value;
          ForLoop<InitSub,0,size-1>::apply(order,points1D,points);
        }
        template <int i>
        struct InitSub
        {
          typedef typename GenericGeometry::SubTopology<Topology,codim,i>::type SubTopology;
          static void apply(const unsigned int order,
                            const std::vector<Field> &points1D,
                            LagrangePoints &points)
          {
            Setup<Topology>::template Init<SubTopology>::template apply<i>(order,points1D,points);
          }
        };
      };
      template <class SubTopology>
      struct Init
      {
        static const unsigned int codimension = dimension - SubTopology::dimension;
        typedef GenericGeometry::ReferenceMappings< Field, dimension > RefMappings;
        typedef typename RefMappings::Container RefMappingsContainer;
        typedef typename RefMappingsContainer::template Codim< codimension >::Mapping Mapping;
        typedef LobattoInnerPoints<Field,SubTopology> InnerPoints;
        template <unsigned int subEntity>
        static void apply(const unsigned int order,
                          const std::vector<Field> &points1D,
                          LagrangePoints &points)
        {
          unsigned int oldSize = points.size();
          unsigned int size = InnerPoints::size(order);
          if (size==0)
            return;
          points.resize(oldSize+size);
          std::vector< LagrangePoint<Field,dimension-codimension> > subPoints(size);

          /*
             std::cout << Topology::name() << " " << SubTopology::name() << " : "
                    << " ( " << codimension <<  " , " << subEntity << " ) "
                    << oldSize << " " << size << std::endl;
           */

          InnerPoints::template setup<dimension-codimension>( points1D,&(subPoints[0]) );

          const RefMappingsContainer &refMappings = RefMappings::container( Topology::id );
          const Mapping &mapping = refMappings.template mapping< codimension >( subEntity );

          LagrangePoint<Field,dimension> *p = &(points.points_[oldSize]);
          for ( unsigned int nr = 0; nr<size; ++nr, ++p)
          {
            /*
               std::cout << "   " << nr << " : "
                      << subPoints[nr].point_ << " -> "
                      << mapping.global( subPoints[nr].point_ )
                      << std::endl;
             */
            p->point_ = mapping.global( subPoints[nr].point_ );
            p->localKey_ = LocalKey( subEntity, codimension, nr );
            #ifndef NDEBUG
            bool test = GenericGeometry::ReferenceElement<Topology,Field>::checkInside(p->point_);
            if (!test)
              std::cout << "not inside" << std::endl;
            // assert( test );
            #endif
          }
        }
      };
    };

    template< class Topology >
    static const LagrangePoints &lagrangePoints ( const Key &order )
    {
      LagrangePoints *lagrangePoints = new LagrangePoints( order, 0 );
      LobattoPoints<Field> points1D(order);
      ForLoop<Setup<Topology>::template InitCodim,0,dimension>::apply(order,points1D.points_,*lagrangePoints);
      return *lagrangePoints;
    }

    template< class T >
    static const LocalCoefficients &localCoefficients ( const Key &order )
    {
      return lagrangePoints< T >( order );
    }

    static void release ( const LagrangePoints &lagrangePoints )
    {
      delete &lagrangePoints;
    }

    template< class Topology >
    static bool supports ( const Key &order )
    {
      const bool isSimplex = GenericGeometry::IsSimplex< Topology >::value;
      const bool isGeneralizedPrism = GenericGeometry::IsGeneralizedPrism< Topology >::value;
      return ( isSimplex || isGeneralizedPrism);
    }
  };

  //
  // LobattoBasisProvider
  // ---------------------

  template< int dim, class SF, class CF = typename ComputeField< SF, 512 >::Type >
  struct LobattoBasisProvider
    : public BasisProvider<
          LagrangeBasisCreator< dim, LobattoPointsCreator, SF, CF > >
  {};
}
#endif // DUNE_LOBATTOBASIS_HH
