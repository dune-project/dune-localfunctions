# $Id: Makefile.am 4150 2008-05-15 16:12:46Z christi $

# which tests to run
TESTS = testfem lagrangeshapefunctiontest

# programs just to build when "make check" is used
check_PROGRAMS = $(TESTS)

# output code coverage
#AM_CXXFLAGS = @AM_CPPFLAGS@ -fprofile-arcs -ftest-coverage

# define the programs

testfem_SOURCES = testfem.cc
lagrangeshapefunctiontest_SOURCES = lagrangeshapefunctiontest.cc

include $(top_srcdir)/am/global-rules
