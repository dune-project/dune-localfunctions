// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_RAVIARTTHOMASBASIS_HH
#define DUNE_RAVIARTTHOMASBASIS_HH
#include <fstream>
#include <utility>

#include "generic.hh"
#include "raviartthomascompute.hh"

namespace Dune
{
  template< int dim, class SF, class CF >
  struct RaviartThomasBasisCreator
  {
    typedef unsigned int Key;
    typedef typename GenericGeometry::SimplexTopology< dim >::type SimplexTopology;

    typedef VirtualMonomialBasis<dim,SF> MBasis;
    typedef SF StorageField;
    typedef AlgLib::MultiPrecision< Precision<CF>::value > ComputeField;
    static const int dimension = dim;
    typedef PolynomialBasisWithMatrix<StandardEvaluator<MBasis>,SparseCoeffMatrix<StorageField,dim> > Basis;
    typedef LocalCoefficientsContainer LocalCoefficients;

    #define RTLAGRANGE 0
#if RTLAGRANGE
    // typedef LagrangePointsCreator< ComputeField, dimension > PointsSetCreator;
    typedef LobattoPointsCreator< ComputeField, dimension > PointsSetCreator;
    typedef RaviartThomasLagrangeInterpolation< ComputeField, PointsSetCreator > LocalInterpolationImpl;
#else
    typedef RaviartThomasL2Interpolation< ComputeField, dimension > LocalInterpolationImpl;
#endif
    typedef RaviartThomasInterpolation< ComputeField, dimension, LocalInterpolationImpl > LocalInterpolation;

    template< class Topology >
    static const LocalInterpolationImpl *localInterpolationImpl ( const Key &order )
    {
      FieldMatrix<ComputeField,dimension+1,dimension> normal;
      for (int f=0; f<dimension+1; ++f)
        normal[f] = GenericGeometry::ReferenceElement<Topology,ComputeField>::integrationOuterNormal(f);
#if RTLAGRANGE
      const typename PointsSetCreator::LagrangePoints &points = PointsSetCreator::template lagrangePoints< Topology >( order+dimension );
      LocalInterpolationImpl *interpolation = new LocalInterpolationImpl(order,points,normal);
#else
      LocalInterpolationImpl *interpolation = new LocalInterpolationImpl(order,normal);
#endif
      return interpolation;
    }
    template< class Topology >
    static const LocalInterpolation &localInterpolation ( const Key &order )
    {
      return * (new LocalInterpolation( localInterpolationImpl<Topology>(order) )) ;
    }
    template <class Topology>
    static Basis &basis(unsigned int order)
    {
      typedef RTPreMatrix<Topology,ComputeField> PreMatrix;
      const MBasis &_mBasis = MonomialBasisProvider<dimension,StorageField>::template basis<Topology>(order+1);
      Basis *basis = new Basis(_mBasis);
      const LocalInterpolationImpl *interpolation = localInterpolationImpl<Topology>(order);
      LocalInterpolation interpolWrapper(interpolation);
      RaviartThomasMatrix<Topology,PreMatrix,LocalInterpolation> matrix(interpolWrapper);
      basis->fill(matrix);
      {
        typedef MultiIndex< dimension > MIField;
        typedef VirtualMonomialBasis<dim,MIField> MBasisMI;
        typedef PolynomialBasisWithMatrix<StandardEvaluator<MBasisMI>,SparseCoeffMatrix<StorageField,dimension> > BasisMI;
        const MBasisMI &_mBasisMI = MonomialBasisProvider<dimension,MIField>::template basis<Topology>(order+1);
        BasisMI basisMI(_mBasisMI);
        basisMI.fill(matrix);
        std::stringstream name;
        name << "rt_" << Topology::name() << "_p" << order;
        std::ofstream out(name.str().c_str());
        basisPrint<0>(out,basisMI);
        /*
           const LocalCoefficients &keys = localCoefficients<Topology>(order);
           for (int index=0;index<keys.size();++index)
           std::cout << index << " -> "
                    << " (codim = " << keys.localKey(index).codim() << ", "
                    << "subentity = " << keys.localKey(index).subEntity() << ", "
                    << "index = " << keys.localKey(index).index() << "):" << std::endl;
         */
      }
      return *basis;
    }
    template< class Topology >
    static const LocalCoefficients &localCoefficients ( const Key &order )
    {
      FieldMatrix<ComputeField,dimension+1,dimension> normal;
      for (int f=0; f<dimension+1; ++f)
        normal[f] = GenericGeometry::ReferenceElement<Topology,ComputeField>::integrationOuterNormal(f);
      const LocalInterpolationImpl *interpolation = localInterpolationImpl<Topology>(order);
      LocalCoefficientsContainer *localKeys = new LocalCoefficientsContainer(*interpolation);
      delete interpolation;
      return *localKeys;
    }

    static void release ( const LocalInterpolation &localInterpolation )
    {
      delete &localInterpolation;
    }
    static void release ( const Basis &basis )
    {
      delete &basis;
    }
    static void release ( const LocalCoefficients &localCoefficients )
    {
      delete &localCoefficients;
    }
    template< class Topology >
    static bool supports ( const Key &key )
    {
      return GenericGeometry::IsSimplex<Topology>::value;
    }
  };

  template< int dim, class SF, class CF = typename ComputeField< SF, 512 >::Type >
  struct RaviartThomasBasisProvider
    : public BasisProvider<RaviartThomasBasisCreator<dim,SF,CF> >
  {};

}
#endif // DUNE_RAVIARTTHOMASBASIS_HH
