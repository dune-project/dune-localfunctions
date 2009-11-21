orthonormalbasisdir = $(includedir)/dune/finiteelements/generic/orthonormalbasis
orthonormalbasis_HEADERS = orthonormalbasis.hh orthonormalcompute.hh
check_PROGRAMS = test
EXTRA_PROGRAMS = test-all

LIBS = $(DUNE_LDFLAGS) $(DUNE_LIBS)

TOPOLOGY = Prism< Pyramid< Pyramid< Point > > >
SETTINGS = -DTOPOLOGY="$(TOPOLOGY)"

test_SOURCES = test.cc
test_CPPFLAGS = $(DUNE_CPPFLAGS) $(AM_CPPFLAGS) $(SETTINGS)
test_all_SOURCES = test.cc
test_all_CPPFLAGS = $(DUNE_CPPFLAGS) $(AM_CPPFLAGS)
LDADD = $(LIBS)

include $(top_srcdir)/am/global-rules


