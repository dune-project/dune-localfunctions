// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LAGRANGEBASIS_HH
#define DUNE_LAGRANGEBASIS_HH

#include <fstream>
#include <dune/alglib/multiprecision.hh>
#include <dune/alglib/matrix.hh>

#include <dune/finiteelements/lagrangebasis/interpolation.hh>
#include <dune/finiteelements/multiindex.hh>
#include <dune/finiteelements/generic/monomialbasis.hh>
#include <dune/finiteelements/generic/basisprovider.hh>
#include <dune/finiteelements/basisprint.hh>
#include <dune/finiteelements/generic/polynomialbasis.hh>

namespace Dune
{

  // LagrangeMatrix
  // --------------

  template< class Topology, class scalar_t,
      class LPCreator >
  struct LagrangeMatrix
  {
    static const unsigned int dimension = Topology::dimension;

    typedef LPCreator LagrangePointsCreator;
    typedef Dune::AlgLib::Matrix< scalar_t > mat_t;
    // typedef Dune::LagrangePointsCreator< scalar_t, dimension > LagrangePointsCreator;
    typedef LocalLagrangeInterpolationCreator< LagrangePointsCreator > LocalInterpolationCreator;
    typedef typename LocalInterpolationCreator::LocalInterpolation LocalInterpolation;

    explicit LagrangeMatrix( const unsigned int order )
    {
      Dune::MonomialBasis< Topology, scalar_t > basis( order );
      const LocalInterpolation &localInterpolation
        = LocalInterpolationCreator::template localInterpolation< Topology >( order );
      localInterpolation.interpolate( basis, matrix_ );
      LocalInterpolationCreator::release( localInterpolation );
      matrix_.invert();
    }

    unsigned int colSize( const unsigned int row ) const
    {
      return matrix_.cols();
    }

    unsigned int rowSize () const
    {
      return matrix_.rows();
    }

    Dune::FieldMatrix< scalar_t, 1, 1 > operator() ( int r, int c ) const
    {
      return matrix_( c, r );
    }

    void print ( std::ostream &out, const unsigned int N = rowSize() ) const
    {
      for( unsigned int i = 0; i < N; ++i )
      {
        out << "Polynomial : " << i << std::endl;
        for( unsigned int j = 0; j <colSize( i ); ++j )
        {
          double v = matrix_( j, i ).toDouble();
          if( fabs( v ) < 1e-20 )
            out << 0 << "\t\t" << std::flush;
          else
          {
            Dune::AlgLib::MultiPrecision< 128 > v = matrix_( j, i );
            out << v << "\t\t" << std::flush;
          }
        }
        out << std::endl;
      }
    }

  private:
    mat_t matrix_;
  };



  // LagrangeBasisCreator
  // --------------------

  template< unsigned int dim,
      template <class Field,unsigned int> class LPCreator,
      class SF, class CF = typename ComputeField< SF, 512 >::Type >
  struct LagrangeBasisCreator
  {
    static const unsigned int dimension = dim;

    typedef SF StorageField;

    typedef StorageField BasisField;

    typedef Dune::MonomialBasisProvider< dimension, BasisField > MonomialBasisProvider;
    typedef typename MonomialBasisProvider::Basis MonomialBasis;

    typedef unsigned int Key;

    typedef StandardEvaluator< MonomialBasis > Evaluator;
    typedef AlgLib::MultiPrecision< Precision< CF >::value > ComputeField;
    typedef PolynomialBasisWithMatrix< Evaluator, SparseCoeffMatrix< StorageField, 1 > > Basis;

    typedef LPCreator<ComputeField,dim> LPointsCreator;

    template< class Topology >
    static const Basis &basis ( const Key &order )
    {
      const MonomialBasis &monomialBasis = MonomialBasisProvider::template basis< Topology >( order );
      Basis *basis = new Basis( monomialBasis );
      LagrangeMatrix< Topology, ComputeField, LPointsCreator> matrix( order );
      basis->fill( matrix );
      {
        typedef MultiIndex< dimension > MIField;
        typedef VirtualMonomialBasis<dim,MIField> MBasisMI;
        typedef PolynomialBasisWithMatrix<StandardEvaluator<MBasisMI>,SparseCoeffMatrix<StorageField,dimension> > BasisMI;
        const MBasisMI &_mBasisMI = Dune::MonomialBasisProvider<dimension,MIField>::template basis<Topology>(order);
        BasisMI basisMI(_mBasisMI);
        basisMI.fill(matrix);
        std::stringstream name;
        name << "rt_" << Topology::name() << "_p" << order;
        std::ofstream out(name.str().c_str());
        basisPrint<0>(out,basisMI);
      }
      return *basis;
    }

    static void release ( const Basis &basis )
    {
      delete &basis;
    }
  };



}

#endif // #ifndef DUNE_ORTHONORMALBASIS_HH
