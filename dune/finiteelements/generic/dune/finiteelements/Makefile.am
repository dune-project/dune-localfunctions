finiteelementsdir = $(includedir)/dune/finiteelements
finiteelements_HEADERS = lagrangefiniteelement.hh

SUBDIRS = generic orthonormalbasis lagrangebasis raviartthomas common quadrature
AM_CPPFLAGS = $(DUNE_CPPFLAGS) -I$(top_srcdir)/alglib

include $(top_srcdir)/am/global-rules
