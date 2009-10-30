lagrangebasisdir = $(includedir)/dune/finiteelements/generic/lagrangebasis
lagrangebasis_HEADERS = emptypoints.hh        lagrangebasis.hh         \
                        equidistantpoints.hh  lagrangecoefficients.hh  \
                        interpolation.hh      lobattopoints.hh

check_PROGRAMS = test test-all


TOPOLOGY = Pyramid< Pyramid< Pyramid< Point > > >
SETTINGS = -DTOPOLOGY="$(TOPOLOGY)"

test_SOURCES = test.cc
test_CPPFLAGS = $(DUNE_CPPFLAGS) $(GMP_CPPFLAGS) $(AM_CPPFLAGS) $(SETTINGS)
test_all_SOURCES = test.cc
test_all_CPPFLAGS = $(DUNE_CPPFLAGS) $(GMP_CPPFLAGS) $(AM_CPPFLAGS)

LDADD = $(DUNE_LDFLAGS) $(DUNE_LIBS) $(GMP_LIBS)

include $(top_srcdir)/am/global-rules


