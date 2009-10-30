AC_DEFUN([DUNE_PATH_ALGLIB],[
  AC_REQUIRE([AC_PROG_CXX])
  AC_REQUIRE([DUNE_PATH_GMP])

  HAVE_ALGLIB=no

  AC_ARG_WITH(alglib,
    AS_HELP_STRING([--with-alglib=PATH],[directory to alglib4dune]))

  ac_save_PKG_CONFIG_PATH="$PKG_CONFIG_PATH"
  AS_IF([test x$with_alglib != x],[PKG_CONFIG_PATH="$with_alglib:$PKG_CONFIG_PATH"])
  AS_IF([pkg-config --atleast-version=1.0 alglib],[
    HAVE_ALGLIB="yes (Version `pkg-config --modversion alglib`)"
    ALGLIB_CPPFLAGS="`pkg-config --cflags alglib` -DENABLE_ALGLIB"
    ALGLIB_LIBS="`pkg-config --libs alglib`"
  ])
  PKG_CONFIG_PATH="$ac_save_PKG_CONFIG_PATH"

  AS_IF([test "$HAVE_ALGLIB" != no],[
    AC_LANG_PUSH([C++])
    ac_save_CPPFLAGS="$CPPFLAGS"
    CPPFLAGS="$CPPFLAGS $ALGLIB_CPPFLAGS"
    AC_CHECK_HEADER([alglib/amp.h],[],[HAVE_ALGLIB=no])
    CPPFLAGS="$ac_save_CPPFLAGS"
    AC_LANG_POP
  ])

  AS_IF([test "$HAVE_ALGLIB" != no],[
    AC_DEFINE([HAVE_ALGLIB],[ENABLE_ALGLIB],[Was AlgLib found and ALGLIB_CPPFLAGS used?])
    AC_SUBST([ALGLIB_CPPFLAGS],[$ALGLIB_CPPFLAGS])
    AC_SUBST([ALGLIB_LIBS],[$ALGLIB_LIBS])
  ])

  AM_CONDITIONAL(ALGLIB,[test "$HAVE_ALGLIB" != no])
  DUNE_ADD_SUMMARY_ENTRY([AlgLib],[$HAVE_ALGLIB])
])
