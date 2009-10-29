raviartthomasbasisdir = $(includedir)/dune/finiteelements/generic/raviartthomasbasis
raviartthomasbasis_HEADERS = interpolation.hh  raviartthomasbasis.hh  raviartthomasprebasis.hh
check_PROGRAMS = test test-all

TOPOLOGY = Pyramid< Pyramid< Point > >
SETTINGS = -DTOPOLOGY="$(TOPOLOGY)"

test_SOURCES = test.cc
test_CPPFLAGS = $(DUNE_CPPFLAGS) $(AM_CPPFLAGS) $(SETTINGS)
test_all_SOURCES = test.cc
test_all_CPPFLAGS = $(DUNE_CPPFLAGS) $(AM_CPPFLAGS)

LDADD = $(DUNE_LDFLAGS) $(DUNE_LIBS)

include $(top_srcdir)/am/global-rules


