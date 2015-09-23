AC_DEFUN([DUNE_LOCALFUNCTIONS_CHECKS],[
  AM_CONDITIONAL([ALUGRID], [test x$HAVE_ALUGRID = x1])
])

AC_DEFUN([DUNE_LOCALFUNCTIONS_CHECK_MODULE],[
  AC_MSG_NOTICE([Searching for dune-localfunctions...])
  DUNE_CHECK_MODULES([dune-localfunctions], [localfunctions/common/localkey.hh])
])
