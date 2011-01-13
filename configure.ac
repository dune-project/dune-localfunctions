# -*- Autoconf -*-
# $Id: configure.ac 5257 2008-08-01 23:19:24Z mblatt $
# Process this file with autoconf to produce a configure script.
AC_PREREQ(2.50)
DUNE_AC_INIT # gets module version from dune.module file
AC_CONFIG_SRCDIR([dune/localfunctions/common/localbasis.hh])
AM_CONFIG_HEADER([config.h])

# add configure flags needed to create log files for dune-autobuild
DUNE_AUTOBUILD_FLAGS
# we need no more than the standard DE-stuff
# this module depends on dune-common
# this implies checking for [dune-common]
DUNE_CHECK_ALL

# preset variable to path such that #include <dune/...> works
# implicitly set the Dune-flags everywhere
AC_SUBST([DUNE_LOCALFUNCTIONS_ROOT], '$(top_builddir)')
AC_SUBST([LOCAL_LIBS], '$(top_builddir)/lib/liblocalfunctions.la')
AC_SUBST(AM_CPPFLAGS, '$(DUNE_CPPFLAGS) -I$(top_srcdir)')
AC_SUBST(AM_LDFLAGS, $DUNE_LDFLAGS)
LIBS="$DUNE_LIBS"

AC_CONFIG_FILES([
		Makefile
                doc/doxygen/Doxyfile
                doc/doxygen/Makefile
                doc/Makefile
                dune-localfunctions.pc
                dune/localfunctions/common/Makefile
                dune/localfunctions/hierarchical/hierarchicalp2/Makefile
                dune/localfunctions/hierarchical/hierarchicalp2withelementbubble/Makefile
                dune/localfunctions/hierarchical/hierarchicalprismp2/Makefile
                dune/localfunctions/hierarchical/Makefile
                dune/localfunctions/lagrange/Makefile
                dune/localfunctions/lagrange/p0/Makefile
                dune/localfunctions/lagrange/p1/Makefile
                dune/localfunctions/lagrange/p23d/Makefile
                dune/localfunctions/lagrange/pk2d/Makefile
                dune/localfunctions/lagrange/pk3d/Makefile
                dune/localfunctions/lagrange/prismp1/Makefile
                dune/localfunctions/lagrange/prismp2/Makefile
                dune/localfunctions/lagrange/pyramidp1/Makefile
                dune/localfunctions/lagrange/q1/Makefile
                dune/localfunctions/lagrange/q22d/Makefile
                dune/localfunctions/Makefile
                dune/localfunctions/mimetic/Makefile
                dune/localfunctions/monom/Makefile
                dune/localfunctions/orthonormal/Makefile
                dune/localfunctions/rannacherturek/Makefile
                dune/localfunctions/rannacherturek/rannacherturek2d/Makefile
                dune/localfunctions/raviartthomas/Makefile
                dune/localfunctions/raviartthomas/raviartthomas02d/Makefile
                dune/localfunctions/raviartthomas/raviartthomas0q2d/Makefile
                dune/localfunctions/raviartthomas/raviartthomas0q3d/Makefile
                dune/localfunctions/refined/common/Makefile
                dune/localfunctions/refined/Makefile
                dune/localfunctions/refined/refinedp0/Makefile
                dune/localfunctions/refined/refinedp1/Makefile
                dune/localfunctions/test/Makefile
                dune/localfunctions/utility/Makefile
                dune/localfunctions/whitney/edges0.5/Makefile
                dune/localfunctions/whitney/Makefile
                dune/Makefile
                m4/Makefile
		])
AC_OUTPUT
# finally print the summary information
DUNE_SUMMARY_ALL
