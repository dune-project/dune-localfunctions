ravtombasisdir = $(includedir)/dune/finiteelements/vectest
ravtombasis_HEADERS = interpolation.hh  raviartthomasbasis.hh  raviartthomasprebasis.hh
check_PROGRAMS = test test-all

TOPOLOGY = Pyramid< Pyramid< Point > >
SETTINGS = -DTOPOLOGY="$(TOPOLOGY)" -DHAVE_GMP=0 -DHAVE_ALGLIB=1

AM_CPPFLAGS = -I$(top_srcdir)/alglib
test_SOURCES = test.cc
test_CPPFLAGS = $(DUNE_CPPFLAGS) $(AM_CPPFLAGS) $(SETTINGS)
test_all_SOURCES = test.cc
test_all_CPPFLAGS = $(DUNE_CPPFLAGS) $(AM_CPPFLAGS)

LDADD = $(DUNE_LDFLAGS) $(DUNE_LIBS) -lgmpxx -lgmp -lmpfr

include $(top_srcdir)/am/global-rules


