# $Id: Makefile.am 4150 2008-05-15 16:12:46Z christi $

if DUNE_GRID
if ALUGRID
  ALUPROG = refinedp1test
endif
endif

# which tests to run
TESTS = lagrangeshapefunctiontest		\
	monomshapefunctiontest			\
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
testgenericfem_CPPFLAGS = $(AM_CPPFLAGS) $(GMP_CPPFLAGS) $(ALGLIB_CPPFLAGS)
testgenericfem_LDADD = $(GMP_LIBS) $(GMP_LDFLAGS) $(ALGLIB_LIBS) $(ALGLIB_LDFLAGS)
lagrangeshapefunctiontest_SOURCES = lagrangeshapefunctiontest.cc
monomshapefunctiontest_SOURCES = monomshapefunctiontest.cc

virtualshapefunctiontest_SOURCES = virtualshapefunctiontest.cc
virtualshapefunctiontest_CPPFLAGS = $(AM_CPPFLAGS)

if DUNE_GRID
if ALUGRID
refinedp1test_SOURCES = refinedp1test.cc
refinedp1test_CPPFLAGS = $(AM_CPPFLAGS)		\
	$(DUNEMPICPPFLAGS)			\
	$(ALUGRID_CPPFLAGS)
refinedp1test_LDFLAGS = $(AM_LDFLAGS)		\
	$(DUNEMPILDFLAGS)			\
	$(ALUGRID_LDFLAGS)
refinedp1test_LDADD =				\
	$(ALUGRID_LDFLAGS) $(ALUGRID_LIBS)	\
	$(DUNEMPILIBS)				\
	$(LDADD)
endif
endif

TESTS += test-edges0.5
test_edges0_5_SOURCES = test-edges0.5.cc

TESTS += test-localfe
test_localfe_SOURCES = test-localfe.cc

TESTS += test-monom
test_monom_SOURCES = test-monom.cc

TESTS += test-pk2d
test_pk2d_SOURCES = test-pk2d.cc

TESTS += test-power-monom
test_power_monom_SOURCES = test-power-monom.cc

TESTS += test-q1
test_q1_SOURCES = test-q1.cc

TESTS += test-q22d
test_q22d_SOURCES = test-q22d.cc

TESTS += test-q2
test_q2_SOURCES = test-q2.cc

include $(top_srcdir)/am/global-rules
