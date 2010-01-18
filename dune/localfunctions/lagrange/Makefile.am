lagrangedir = $(includedir)/dune/localfunctions/lagrange
lagrange_HEADERS = emptypoints.hh        lagrangebasis.hh         \
                   equidistantpoints.hh  lagrangecoefficients.hh  \
                   interpolation.hh      lobattopoints.hh

check_PROGRAMS = test-all
EXTRA_PROGRMAS = test


TOPOLOGY = Pyramid< Pyramid< Pyramid< Point > > >
SETTINGS = -DTOPOLOGY="$(TOPOLOGY)"

test_SOURCES = test.cc
test_CPPFLAGS = $(DUNE_CPPFLAGS) $(GMP_CPPFLAGS) $(ALGLIB_CPPFLAGS) $(AM_CPPFLAGS) $(SETTINGS)
test_all_SOURCES = test.cc
test_all_CPPFLAGS = $(DUNE_CPPFLAGS) $(GMP_CPPFLAGS) $(ALGLIB_CPPFLAGS) $(AM_CPPFLAGS)

LDADD = $(DUNE_LDFLAGS) $(DUNE_LIBS) $(GMP_LIBS) $(ALGLIB_LIBS)

include $(top_srcdir)/am/global-rules
