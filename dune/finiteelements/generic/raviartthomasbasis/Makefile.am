raviartthomasbasisdir = $(includedir)/dune/finiteelements/generic/raviartthomasbasis
raviartthomasbasis_HEADERS = interpolation.hh  raviartthomasbasis.hh  raviartthomasprebasis.hh
check_PROGRAMS = test
EXTRA_PROGRMAS = test-all

TOPOLOGY = Pyramid< Pyramid< Point > >
SETTINGS = -DTOPOLOGY="$(TOPOLOGY)"

test_SOURCES = test.cc
test_CPPFLAGS = $(DUNE_CPPFLAGS) $(AM_CPPFLAGS) $(SETTINGS) $(GMP_CPPFLAGS) $(ALGLIB_CPPFLAGS)
test_all_SOURCES = test.cc
test_all_CPPFLAGS = $(DUNE_CPPFLAGS) $(AM_CPPFLAGS) $(GMP_CPPFLAGS) $(ALGLIB_CPPFLAGS)


LDADD = $(DUNE_LDFLAGS) $(DUNE_LIBS) $(GMP_LIBS) $(ALGLIB_LIBS)

include $(top_srcdir)/am/global-rules


