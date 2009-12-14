orthonormalbasisdir = $(includedir)/dune/localfunctions/generic/orthonormalbasis
orthonormalbasis_HEADERS = orthonormalbasis.hh orthonormalcompute.hh
check_PROGRAMS = test
EXTRA_PROGRAMS = test-all

TOPOLOGY = Prism< Pyramid< Pyramid< Point > > >
SETTINGS = -DTOPOLOGY="$(TOPOLOGY)"

test_SOURCES = test.cc
test_CPPFLAGS = $(DUNE_CPPFLAGS) $(GMP_CPPFLAGS) $(ALGLIB_CPPFLAGS) $(AM_CPPFLAGS) $(SETTINGS)
test_all_SOURCES = test.cc
test_all_CPPFLAGS = $(DUNE_CPPFLAGS) $(GMP_CPPFLAGS) $(ALGLIB_CPPFLAGS) $(AM_CPPFLAGS)

LDADD = $(DUNE_LDFLAGS) $(DUNE_LIBS) $(GMP_LIBS) $(ALGLIB_LIBS)

include $(top_srcdir)/am/global-rules


