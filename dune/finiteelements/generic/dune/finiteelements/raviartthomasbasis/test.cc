// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <dune/finiteelements/common/field.hh>
#include <dune/finiteelements/raviartthomasbasis/raviartthomasbasis.hh>
#include <dune/finiteelements/generic/basisprint.hh>

#if HAVE_ALGLIB
typedef amp::ampf< 128 > StorageField;
typedef amp::ampf< 512 > ComputeField;
#else
#if HAVE_GMP
typedef Dune::GMPField< 128 > StorageField;
typedef Dune::GMPField< 512 > ComputeField;
#else
typedef double StorageField;
typedef double ComputeField;
#endif
#endif

template <class Topology>
bool test(unsigned int order)
{
  bool ret = true;

  for (unsigned int o=order; o>=0; --o)
  {
    std::cout << "Testing " << Topology::name() << " in dimension " << Topology::dimension << " with order " << o << std::endl;
    typedef Dune::RaviartThomasBasisFactory<Topology::dimension,StorageField,ComputeField> BasisFactory;
    const typename BasisFactory::Object &basis = *BasisFactory::template create<Topology>(o);
    std::stringstream name;
    name << "rt_" << Topology::name() << "_p" << o << ".basis";
    std::ofstream out(name.str().c_str());
    Dune::basisPrint<0,BasisFactory>(out,basis);
    BasisFactory::release(&basis);
  }
  if (!ret) {
    std::cout << "   FAILED !" << std::endl;
  }
  return ret;
}
int main ( int argc, char **argv )
{
  using namespace Dune;
  using namespace GenericGeometry;

  if( argc < 2 )
  {
    std::cerr << "Usage: " << argv[ 0 ] << " <p>" << std::endl;
    return 2;
  }

  const unsigned int order = atoi( argv[ 1 ] );
#ifdef TOPOLOGY
  return (test<TOPOLOGY>(order) ? 0 : 1 );
#else
  bool tests = true;
  tests &= test<Pyramid<Point> > (order);

  tests &= test<Pyramid<Pyramid<Point> > >(order);

  tests &= test<Pyramid<Pyramid<Pyramid<Point> > > >(order);

  tests &= test<Pyramid<Pyramid<Pyramid<Pyramid<Point> > > > >(order);
  return (tests ? 0 : 1);
#endif // TOPOLOGY
}
