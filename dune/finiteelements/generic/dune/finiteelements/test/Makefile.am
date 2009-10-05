TOPOLOGY = Pyramid< Pyramid< Point > >
# TOPOLOGY = Prism< Prism< Point > >

check_PROGRAMS = test test-lagrange test-jacobian test-monoms test-onb test-vec test-dglagrange test-lobatto test-quadrature

SETTINGS = -DTOPOLOGY="$(TOPOLOGY)"

AM_CPPFLAGS = $(DUNE_CPPFLAGS) -I$(top_srcdir)/alglib
LDADD = $(DUNE_LDFLAGS) $(DUNE_LIBS) -lmpfr -lgmp

test_SOURCES = test.cc
test_CPPFLAGS = $(AM_CPPFLAGS) $(SETTINGS)

test_vec_SOURCES = test-vec.cc
test_vec_CPPFLAGS = $(AM_CPPFLAGS) $(SETTINGS)

test_quadrature_SOURCES = test-quadrature.cc
test_quadrature_CPPFLAGS = $(AM_CPPFLAGS) $(SETTINGS)

test_jacobian_SOURCES = test-jacobian.cc
test_jacobian_CPPFLAGS = $(AM_CPPFLAGS) $(SETTINGS)

test_monoms_SOURCES = test-monoms.cc
test_monoms_CPPFLAGS = $(AM_CPPFLAGS) $(SETTINGS)

test_lagrange_SOURCES = test-space.cc
test_lagrange_CPPFLAGS = -DLAGRANGE $(AM_CPPFLAGS) $(ALL_PKG_CPPFLAGS)
test_lagrange_LDADD = $(LDADD) $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS)

test_lobatto_SOURCES = test-space.cc
test_lobatto_CPPFLAGS = -DLOBATTO $(AM_CPPFLAGS) $(ALL_PKG_CPPFLAGS)
test_lobatto_LDADD = $(LDADD) $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS)

test_dglagrange_SOURCES = test-space.cc
test_dglagrange_CPPFLAGS = -DDGLAGRANGE $(AM_CPPFLAGS) $(ALL_PKG_CPPFLAGS)
test_dglagrange_LDADD = $(LDADD) $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS)

test_onb_SOURCES = test-space.cc
test_onb_CPPFLAGS = -DONB $(AM_CPPFLAGS) $(ALL_PKG_CPPFLAGS)
test_onb_LDADD = $(LDADD) $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS)

include $(top_srcdir)/am/global-rules
