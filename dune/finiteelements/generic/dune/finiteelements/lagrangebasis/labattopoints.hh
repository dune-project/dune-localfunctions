// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
template <class Field>
struct LabattoPoints1D
{
  LabattoPoints1D(unsigned int order) {}
  const Field &point(int i)
  {
    return points_[i];
  }
  std::vector<Field> points_;
};

template <class Field,class Topology>
struct LabattoInnerPoints;

template <class Field>
struct LabattoInnerPoints<Field,Point>
{
  static const unsigned int dim = 0;
  static void setup(const unsigned int order,
                    LagrangePoint< Field, dim > *points )
  {}
};
template <class Field,class Base>
struct LabattoInnerPoints<Field,Pyramid<Base> >
{
  static const unsigned int dim = Base::dimension+1;
  static void setup(const unsigned int order,
                    LagrangePoint< Field, dim > *points )
  {
    // ...
  }
};
template <class Field,class Base>
struct LabattoInnerPoints<Field,Prism<Base> >
{
  static const unsigned int dim = Base::dimension+1;
  static void setup(const unsigned int order,
                    LagrangePoint< Field, dim > *points )
  {}
};
template <class Field>

template <class Field,int dimension>
struct LabattoPoints
{
  template <class Topology>
  struct Creator
  {
    template <unsigned int dim>
    static void apply( unsigned int order,
                       LagrangePoint< Field, dim > *points )
    {
      LabattoInnerPoints<Field,Topology>::setup(order,points);
    }
  };
  static void setup(unsigned int order,
                    unsigned int topoId,
                    LagrangePoint< Field, dimension > *points )
  {
    // get Reference element refElem
    for (int c=dimension; c>=0; --c)
    {
      for (int i=0; i<refElem.size(c); ++i)
      {
        // get subEntity reference element and mapping
        LagrangePoint<Field,dimension-c> pSub;
        GenericGeometry::IfTopology<Creator,dimension>::apply(subTopoId,order,points);
        // map pSub to refElem and add to points
      }
    }
  }
};
