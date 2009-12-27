# $Id: Makefile.am 4150 2008-05-15 16:12:46Z christi $

if ALUGRID
  ALUPROG = refinedp1test
endif

# which tests to run
TESTS = testfem edges02dshapefunctiontest lagrangeshapefunctiontest \
	monomshapefunctiontest			\
	test-edges03dinterpolate		\
	virtualshapefunctiontest		\
	$(ALUPROG)


# programs just to build when "make check" is used
check_PROGRAMS = $(TESTS)

# output code coverage
#AM_CXXFLAGS = @AM_CPPFLAGS@ -fprofile-arcs -ftest-coverage

# define the programs

testfem_SOURCES = testfem.cc
edges02dshapefunctiontest_SOURCES = edges02dshapefunctiontest.cc
lagrangeshapefunctiontest_SOURCES = lagrangeshapefunctiontest.cc
monomshapefunctiontest_SOURCES = monomshapefunctiontest.cc
test_edges03dinterpolate_SOURCES = test-edges03dinterpolate.cc

virtualshapefunctiontest_SOURCES = virtualshapefunctiontest.cc
virtualshapefunctiontest_CPPFLAGS = $(AM_CPPFLAGS)

if ALUGRID
refinedp1test_SOURCES = refinedp1test.cc
refinedp1test_CXXFLAGS = $(ALUGRID_CPPFLAGS) $(MPI_CPPFLAGS) $(AM_CPPFLAGS)
refinedp1test_LDADD = $(ALUGRID_LDFLAGS) $(ALUGRID_LIBS) $(MPI_LDFLAGS)
endif




include $(top_srcdir)/am/global-rules
