orthonormalbasisdir = $(includedir)/dune/finiteelements/orthonormalbasis
orthonormalbasis_HEADERS = orthonormalbasis.hh orthonormalcompute.hh l2interpolation.hh
check_PROGRAMS = example test test-all

LIBS = $(DUNE_LDFLAGS) $(DUNE_LIBS) -lgmp -lmpfr

TOPOLOGY = Prism< Pyramid< Pyramid< Point > > >
SETTINGS = -DTOPOLOGY="$(TOPOLOGY)"

AM_CPPFLAGS = -I$(top_srcdir)/alglib
example_SOURCES = example.cc
example_CPPFLAGS = $(DUNE_CPPFLAGS) $(AM_CPPFLAGS)
test_SOURCES = test.cc
test_CPPFLAGS = $(DUNE_CPPFLAGS) $(AM_CPPFLAGS) $(SETTINGS)
test_all_SOURCES = test.cc
test_all_CPPFLAGS = $(DUNE_CPPFLAGS) $(AM_CPPFLAGS)
LDADD = $(LIBS)

include $(top_srcdir)/am/global-rules


