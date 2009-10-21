finiteelementsdir = $(includedir)/dune/finiteelements
finiteelements_HEADERS = lagrangefiniteelement.hh \
                         orthonormalfiniteelement.hh \
                         raviartthomasfiniteelement.hh
SUBDIRS = generic orthonormalbasis lagrangebasis raviartthomasbasis common quadrature
AM_CPPFLAGS = $(DUNE_CPPFLAGS) -I$(top_srcdir)/alglib

include $(top_srcdir)/am/global-rules
