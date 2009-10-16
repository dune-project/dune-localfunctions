finiteelementsdir = $(includedir)/dune/finiteelements
finiteelements_HEADERS = multiindex.hh dglocalcoefficients.hh

SUBDIRS = test generic orthonormalbasis lagrangebasis raviartthomas global common
AM_CPPFLAGS = $(DUNE_CPPFLAGS) -I$(top_srcdir)/alglib

include $(top_srcdir)/am/global-rules
