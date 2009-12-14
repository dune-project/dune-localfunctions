// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LOBATTOBASIS_HH
#define DUNE_LOBATTOBASIS_HH

#include <fstream>
#include <dune/localfunctions/generic/math/matrix.hh>
#include <dune/localfunctions/generic/math/field.hh>
#include <dune/common/forloop.hh>
#include <dune/localfunctions/generic/common/topologyfactory.hh>
#include <dune/grid/genericgeometry/referenceelements.hh>
#include <dune/grid/genericgeometry/referencemappings.hh>

#include <dune/localfunctions/generic/quadrature/lobattoquadrature.hh>
#include <dune/localfunctions/generic/lagrangebasis/lagrangecoefficients.hh>
#include <dune/localfunctions/generic/lagrangebasis/emptypoints.hh>

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
#if HAVE_ALGLIB
      typedef amp::ampf< Precision< Field >::value > MPField;
#else
      typedef Field MPField;
#endif
      GenericGeometry::LobattoPoints<MPField> lobatto(order+1);

      for (unsigned int i=1; i<order; ++i) {
        points_[i-1] = field_cast<Field>(lobatto.point(i));
      }
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

  template< class F, unsigned int dim >
  struct LobattoPointSet : public EmptyPointSet<F,dim>
  {
    // friend class LagrangeCoefficientsFactory<LobattoPointSet,dim,F>;
    static const unsigned int dimension = dim;
    typedef F Field;
    typedef EmptyPointSet<F,dim> Base;
    typedef typename Base::LagrangePoint Point;
    LobattoPointSet(unsigned int order)
      : Base(order)
    {}

    template< class Topology >
    bool build ( )
    {
      unsigned int order = Base::order();
      LobattoPoints<Field> points1D(order);
      ForLoop<Setup<Topology>::template InitCodim,0,dimension>::
      apply(order,points1D.points_,points_);
      return true;
    }
    template< class Topology >
    static bool supports ( unsigned int order )
    {
      if ( GenericGeometry::IsSimplex< Topology >::value ||
           GenericGeometry::IsGeneralizedPrism< Topology >::value )
        return true;
      else
        return false;
    }
  protected:
    using Base::points_;
  private:
    template <class Topology>
    struct Setup
    {
      template <int pdim>
      struct InitCodim
      {
        static const unsigned int codim = dimension-pdim;
        static void apply(const unsigned int order,
                          const std::vector<Field> &points1D,
                          std::vector<Point> &points)
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
                            std::vector<Point> &points)
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
                          std::vector<Point> &points)
        {
          unsigned int oldSize = points.size();
          unsigned int size = InnerPoints::size(order);
          if (size==0)
            return;
          points.resize(oldSize+size);
          std::vector< LagrangePoint<Field,dimension-codimension> > subPoints(size);

          InnerPoints::template setup<dimension-codimension>( points1D,&(subPoints[0]) );

          const RefMappingsContainer &refMappings = RefMappings::container( Topology::id );
          const Mapping &mapping = refMappings.template mapping< codimension >( subEntity );

          LagrangePoint<Field,dimension> *p = &(points[oldSize]);
          for ( unsigned int nr = 0; nr<size; ++nr, ++p)
          {
            p->point_ = mapping.global( subPoints[nr].point_ );
            p->localKey_ = LocalKey( subEntity, codimension, nr );
            #ifndef NDEBUG
            bool test = GenericGeometry::ReferenceElement<Topology,Field>::checkInside(p->point_);
            if (!test)
              std::cout << "not inside" << std::endl;
            #endif
          }
        }
      };
    };

  };
}
#endif // DUNE_LOBATTOBASIS_HH
