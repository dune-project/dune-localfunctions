# $Id: Makefile.am 4150 2008-05-15 16:12:46Z christi $

# which tests to run
TESTS = testfem edges02dshapefunctiontest lagrangeshapefunctiontest \
	monomshapefunctiontest			\
	test-edges03dinterpolate		\
	virtualshapefunctiontest

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

include $(top_srcdir)/am/global-rules
