# $Id: duneproject 5232 2008-07-03 09:34:11Z christi $

# we need the module file to be able to build via dunecontrol
EXTRA_DIST=dune.module

# don't follow the full GNU-standard
# we need automake 1.5
AUTOMAKE_OPTIONS = foreign 1.5

SUBDIRS = dune doc m4

# pass most important options when "make distcheck" is used
DISTCHECK_CONFIGURE_FLAGS = --with-dune-common="$(DUNE_COMMON_ROOT)" CXX="$(CXX)" CC="$(CC)"

include $(top_srcdir)/am/top-rules
include $(top_srcdir)/am/global-rules
