# which tests to run
TESTS = lagrangeshapefunctiontest		\
	monomialshapefunctiontest			\
	virtualshapefunctiontest		\
	testgenericfem              \
	$(ALUPROG)


# programs just to build when "make check" is used
check_PROGRAMS = $(TESTS)

# output code coverage
#AM_CXXFLAGS = @AM_CPPFLAGS@ -fprofile-arcs -ftest-coverage

noinst_HEADERS =				\
	geometries.hh				\
	test-fe.hh				\
	test-localfe.hh

# define the programs

testgenericfem_SOURCES = testgenericfem.cc
testgenericfem_CPPFLAGS = $(AM_CPPFLAGS) $(GMP_CPPFLAGS)
testgenericfem_LDADD = $(GMP_LIBS) $(GMP_LDFLAGS)

lagrangeshapefunctiontest_SOURCES = lagrangeshapefunctiontest.cc
monomialshapefunctiontest_SOURCES = monomialshapefunctiontest.cc

virtualshapefunctiontest_SOURCES = virtualshapefunctiontest.cc
virtualshapefunctiontest_CPPFLAGS = $(AM_CPPFLAGS)

TESTS += test-edges0.5
test_edges0_5_SOURCES = test-edges0.5.cc

TESTS += test-localfe
test_localfe_SOURCES = test-localfe.cc

TESTS += test-monomial
test_monomial_SOURCES = test-monomial.cc

TESTS += test-pk2d
test_pk2d_SOURCES = test-pk2d.cc

TESTS += test-power-monomial
test_power_monomial_SOURCES = test-power-monomial.cc

TESTS += test-q1
test_q1_SOURCES = test-q1.cc

TESTS += test-q2
test_q2_SOURCES = test-q2.cc

# test-lagrange
TESTS += test-lagrange1
LAGRANGE1_TOPOLOGY = Pyramid< Point >
LAGRANGE1_SETTINGS = -DTOPOLOGY="$(LAGRANGE1_TOPOLOGY)"
test_lagrange1_SOURCES = test-lagrange.cc
test_lagrange1_CPPFLAGS = $(DUNE_CPPFLAGS) $(AM_CPPFLAGS) $(LAGRANGE1_SETTINGS) $(GMP_CPPFLAGS)
test_lagrange1_LDADD = $(DUNE_LDFLAGS) $(DUNE_LIBS) $(GMP_LIBS)

TESTS += test-lagrange2
LAGRANGE2_TOPOLOGY = Pyramid< Pyramid< Point > >
LAGRANGE2_SETTINGS = -DTOPOLOGY="$(LAGRANGE2_TOPOLOGY)"
test_lagrange2_SOURCES = test-lagrange.cc
test_lagrange2_CPPFLAGS = $(DUNE_CPPFLAGS) $(AM_CPPFLAGS) $(LAGRANGE2_SETTINGS) $(GMP_CPPFLAGS)
test_lagrange2_LDADD = $(DUNE_LDFLAGS) $(DUNE_LIBS) $(GMP_LIBS)

TESTS += test-lagrange3
LAGRANGE3_TOPOLOGY = Pyramid< Pyramid< Pyramid< Point > > >
LAGRANGE3_SETTINGS = -DTOPOLOGY="$(LAGRANGE3_TOPOLOGY)"
test_lagrange3_SOURCES = test-lagrange.cc
test_lagrange3_CPPFLAGS = $(DUNE_CPPFLAGS) $(AM_CPPFLAGS) $(LAGRANGE3_SETTINGS) $(GMP_CPPFLAGS)
test_lagrange3_LDADD = $(DUNE_LDFLAGS) $(DUNE_LIBS) $(GMP_LIBS)

# test-orthonormal
TESTS += test-orthonormal1
ORTHONORMAL1_TOPOLOGY = Pyramid< Point >
ORTHONORMAL1_SETTINGS = -DTOPOLOGY="$(ORTHONORMAL1_TOPOLOGY)"
test_orthonormal1_SOURCES = test-orthonormal.cc
test_orthonormal1_CPPFLAGS = $(DUNE_CPPFLAGS) $(AM_CPPFLAGS) $(ORTHONORMAL1_SETTINGS) $(GMP_CPPFLAGS)
test_orthonormal1_LDADD = $(DUNE_LDFLAGS) $(DUNE_LIBS) $(GMP_LIBS)

TESTS += test-orthonormal2
ORTHONORMAL2_TOPOLOGY = Pyramid< Pyramid< Point > >
ORTHONORMAL2_SETTINGS = -DTOPOLOGY="$(ORTHONORMAL2_TOPOLOGY)"
test_orthonormal2_SOURCES = test-orthonormal.cc
test_orthonormal2_CPPFLAGS = $(DUNE_CPPFLAGS) $(AM_CPPFLAGS) $(ORTHONORMAL2_SETTINGS) $(GMP_CPPFLAGS)
test_orthonormal2_LDADD = $(DUNE_LDFLAGS) $(DUNE_LIBS) $(GMP_LIBS)

TESTS += test-orthonormal3
ORTHONORMAL3_TOPOLOGY = Prism< Pyramid< Pyramid< Point > > >
ORTHONORMAL3_SETTINGS = -DTOPOLOGY="$(ORTHONORMAL3_TOPOLOGY)"
test_orthonormal3_SOURCES = test-orthonormal.cc
test_orthonormal3_CPPFLAGS = $(DUNE_CPPFLAGS) $(AM_CPPFLAGS) $(ORTHONORMAL3_SETTINGS) $(GMP_CPPFLAGS)
test_orthonormal3_LDADD = $(DUNE_LDFLAGS) $(DUNE_LIBS) $(GMP_LIBS)

# test-raviartthomassimplex
TESTS += test-raviartthomas1simplex
RT1_TOPOLOGY = Pyramid< Point >
RT1_SETTINGS = -DTOPOLOGY="$(RT1_TOPOLOGY)"
test_raviartthomas1simplex_SOURCES = test-raviartthomassimplex.cc
test_raviartthomas1simplex_CPPFLAGS = $(DUNE_CPPFLAGS) $(AM_CPPFLAGS) $(RT1_SETTINGS) $(GMP_CPPFLAGS)
test_raviartthomas1simplex_LDADD = $(DUNE_LDFLAGS) $(DUNE_LIBS) $(GMP_LIBS)

TESTS += test-raviartthomas2simplex
RT2_TOPOLOGY = Pyramid< Pyramid< Point > >
RT2_SETTINGS = -DTOPOLOGY="$(RT2_TOPOLOGY)"
test_raviartthomas2simplex_SOURCES = test-raviartthomassimplex.cc
test_raviartthomas2simplex_CPPFLAGS = $(DUNE_CPPFLAGS) $(AM_CPPFLAGS) $(RT2_SETTINGS) $(GMP_CPPFLAGS)
test_raviartthomas2simplex_LDADD = $(DUNE_LDFLAGS) $(DUNE_LIBS) $(GMP_LIBS)

TESTS += test-raviartthomas3simplex
RT3_TOPOLOGY = Pyramid< Pyramid< Pyramid< Point > > >
RT3_SETTINGS = -DTOPOLOGY="$(RT3_TOPOLOGY)"
test_raviartthomas3simplex_SOURCES = test-raviartthomassimplex.cc
test_raviartthomas3simplex_CPPFLAGS = $(DUNE_CPPFLAGS) $(AM_CPPFLAGS) $(RT3_SETTINGS) $(GMP_CPPFLAGS)
test_raviartthomas3simplex_LDADD = $(DUNE_LDFLAGS) $(DUNE_LIBS) $(GMP_LIBS)

# extra programs
EXTRA_PROGRAMS = \
	test-lagrange-all \
	test-orthonormal-all \
	test-raviartthomassimplex-all

# these are not included in automatic testing due to long compilation
test_lagrange_all_SOURCES = test-lagrange.cc
test_lagrange_all_CPPFLAGS = $(DUNE_CPPFLAGS) $(AM_CPPFLAGS) $(GMP_CPPFLAGS)
test_lagrange_all_LDADD = $(DUNE_LDFLAGS) $(DUNE_LIBS) $(GMP_LIBS)

test_orthonormal_all_SOURCES = test-orthonormal.cc
test_orthonormal_all_CPPFLAGS = $(DUNE_CPPFLAGS) $(AM_CPPFLAGS) $(GMP_CPPFLAGS)
test_orthonormal_all_LDADD = $(DUNE_LDFLAGS) $(DUNE_LIBS) $(GMP_LIBS)

test_raviartthomassimplex_all_SOURCES = test-raviartthomassimplex.cc
test_raviartthomassimplex_all_CPPFLAGS = $(DUNE_CPPFLAGS) $(AM_CPPFLAGS) $(GMP_CPPFLAGS)
test_raviartthomassimplex_all_LDADD = $(DUNE_LDFLAGS) $(DUNE_LIBS) $(GMP_LIBS)

include $(top_srcdir)/am/global-rules

EXTRA_DIST = CMakeLists.txt
