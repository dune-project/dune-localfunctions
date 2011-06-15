SUBDIRS = p0 p1 p23d pk2d pk3d prismp1 prismp2 pyramidp1 q1 q2 q22d q23d

lagrangedir = $(includedir)/dune/localfunctions/lagrange
lagrange_HEADERS = emptypoints.hh        lagrangebasis.hh         \
                   equidistantpoints.hh  lagrangecoefficients.hh  \
                   interpolation.hh      lobattopoints.hh         \
                   p0.hh p1.hh p23d.hh p2.hh pk2d.hh pk3d.hh pk.hh pq22d.hh \
                   prismp1.hh prismp2.hh pyramidp1.hh q1.hh q2.hh q22d.hh q23d.hh \
                   pqkfactory.hh

check_PROGRAMS = test-all
EXTRA_PROGRAMS = test


TOPOLOGY = Pyramid< Pyramid< Pyramid< Point > > >
SETTINGS = -DTOPOLOGY="$(TOPOLOGY)"

test_SOURCES = test.cc
test_CPPFLAGS = $(DUNE_CPPFLAGS) $(GMP_CPPFLAGS) $(ALGLIB_CPPFLAGS) $(AM_CPPFLAGS) $(SETTINGS)
test_all_SOURCES = test.cc
test_all_CPPFLAGS = $(DUNE_CPPFLAGS) $(GMP_CPPFLAGS) $(ALGLIB_CPPFLAGS) $(AM_CPPFLAGS)

LDADD = $(DUNE_LDFLAGS) $(DUNE_LIBS) $(GMP_LIBS) $(ALGLIB_LIBS)

include $(top_srcdir)/am/global-rules
