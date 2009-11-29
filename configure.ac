# -*- Autoconf -*-
# $Id: configure.ac 5257 2008-08-01 23:19:24Z mblatt $
# Process this file with autoconf to produce a configure script.
AC_PREREQ(2.50)
DUNE_AC_INIT # gets module version from dune.module file
AM_INIT_AUTOMAKE(tar-ustar)
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
                dune/finiteelements/common/Makefile
                dune/finiteelements/edger02d/Makefile
                dune/finiteelements/edges02d/Makefile
                dune/finiteelements/edges03d/Makefile
                dune/finiteelements/Makefile
                dune/finiteelements/monom/Makefile
                dune/finiteelements/mimetic/Makefile
                dune/finiteelements/p0/Makefile
                dune/finiteelements/p1/Makefile
                dune/finiteelements/pk2d/Makefile
                dune/finiteelements/pk3d/Makefile
                dune/finiteelements/prismp1/Makefile
                dune/finiteelements/prismp2/Makefile
                dune/finiteelements/q1/Makefile
                dune/finiteelements/q22d/Makefile
                dune/finiteelements/hierarchicalprismp2/Makefile
                dune/finiteelements/refinedp1/Makefile
                dune/finiteelements/rannacher_turek2d/Makefile
                dune/finiteelements/rt02d/Makefile
                dune/finiteelements/rt0q2d/Makefile
                dune/finiteelements/rt0q3d/Makefile
                dune/finiteelements/test/Makefile
                dune/finiteelements/generic/Makefile
                dune/finiteelements/generic/common/Makefile
                dune/finiteelements/generic/math/Makefile
                dune/finiteelements/generic/quadrature/Makefile
                dune/finiteelements/generic/orthonormalbasis/Makefile
                dune/finiteelements/generic/lagrangebasis/Makefile
                dune/finiteelements/generic/raviartthomasbasis/Makefile
                doc/doxygen/Makefile
                doc/doxygen/Doxyfile
                doc/Makefile
                doc/Makefile.dist
                debug/Makefile
                m4/Makefile
                dune-localfunctions.pc
		])
AC_OUTPUT
# finally print the summary information
DUNE_SUMMARY_ALL
