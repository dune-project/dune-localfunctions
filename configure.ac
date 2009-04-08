# -*- Autoconf -*-
# $Id: configure.ac 5257 2008-08-01 23:19:24Z mblatt $
# Process this file with autoconf to produce a configure script.
AC_PREREQ(2.50)
DUNE_AC_INIT # gets module version from dune.module file
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([dune/finiteelements/common/localbasis.hh])
AM_CONFIG_HEADER([config.h])

# add configure flags needed to create log files for dune-autobuild
DUNE_AUTOBUILD_FLAGS
# we need no more than the standard DE-stuff
# this module depends on dune-common
# this implies checking for [dune-common]
DUNE_CHECK_ALL

# preset variable to path such that #include <dune/...> works
AC_SUBST([DUNE_LOCALFUNCTIONS_ROOT], '$(top_builddir)')
AC_SUBST([AM_CPPFLAGS], '-I$(top_srcdir)')
AC_SUBST([LOCAL_LIBS], '$(top_builddir)/lib/liblocalfunctions.la')

# implicitly set the Dune-flags everywhere
AC_SUBST(AM_CPPFLAGS, $DUNE_CPPFLAGS)
AC_SUBST(AM_LDFLAGS, $DUNE_LDFLAGS)
LIBS="$DUNE_LIBS"

AC_CONFIG_FILES([
		Makefile
                dune/Makefile
                dune/finiteelements/Makefile
                dune/finiteelements/common/Makefile
                dune/finiteelements/test/Makefile
                dune/finiteelements/p0/Makefile
                dune/finiteelements/p11d/Makefile
                dune/finiteelements/p12d/Makefile
                dune/finiteelements/p13d/Makefile
                dune/finiteelements/pk2d/Makefile
                dune/finiteelements/q12d/Makefile
                dune/finiteelements/q13d/Makefile
                dune/finiteelements/q22d/Makefile
                dune/finiteelements/rt02d/Makefile
                dune/finiteelements/edger12d/Makefile
                dune/finiteelements/refinedp1/Makefile
                dune/finiteelements/monom/Makefile
		lib/Makefile
                doc/layout/Makefile
                doc/doxygen/Makefile
                doc/doxygen/Doxyfile
                doc/Makefile
                doc/Makefile.dist
                debug/Makefile
                dune-localfunctions.pc
		])
AC_OUTPUT
# finally print the summary information
DUNE_SUMMARY_ALL
