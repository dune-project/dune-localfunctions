AC_DEFUN([DUNE_LOCALFUNCTIONS_CHECKS],[
  AM_CONDITIONAL([DUNE_GRID], [test x"$with_dune_grid" = xyes])
])

AC_DEFUN([DUNE_LOCALFUNCTIONS_CHECK_MODULE],[
  AC_MSG_NOTICE([Searching for dune-localfunctions...])
  DUNE_CHECK_MODULES([dune-localfunctions], [finiteelements/common/localbasis.hh])
])
