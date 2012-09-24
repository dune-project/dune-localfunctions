SUBDIRS = p0 p1 p23d pk1d pk2d pk3d prismp1 prismp2 pyramidp1 pyramidp2 q1 q2 qk

lagrangedir = $(includedir)/dune/localfunctions/lagrange
lagrange_HEADERS = emptypoints.hh        lagrangebasis.hh         \
                   equidistantpoints.hh  lagrangecoefficients.hh  \
                   interpolation.hh \
                   p0.hh p1.hh p23d.hh p2.hh pk1d.hh pk2d.hh pk3d.hh pk.hh pq22d.hh \
                   prismp1.hh prismp2.hh pyramidp1.hh pyramidp2.hh q1.hh q2.hh \
                   qk.hh \
                   pqkfactory.hh

check_PROGRAMS = test-all
EXTRA_PROGRAMS = test


TOPOLOGY = Pyramid< Pyramid< Pyramid< Point > > >
SETTINGS = -DTOPOLOGY="$(TOPOLOGY)"

test_SOURCES = test.cc
test_CPPFLAGS = $(DUNE_CPPFLAGS) $(GMP_CPPFLAGS) $(AM_CPPFLAGS) $(SETTINGS)
test_all_SOURCES = test.cc
test_all_CPPFLAGS = $(DUNE_CPPFLAGS) $(GMP_CPPFLAGS) $(AM_CPPFLAGS)

LDADD = $(DUNE_LDFLAGS) $(DUNE_LIBS) $(GMP_LIBS)

include $(top_srcdir)/am/global-rules
