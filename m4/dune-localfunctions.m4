AC_DEFUN([DUNE_LOCALFUNCTIONS_CHECKS],[])

AC_DEFUN([DUNE_LOCALFUNCTIONS_CHECK_MODULE],[
  AC_MSG_NOTICE([Searching for dune-localfunctions...])
  DUNE_CHECK_MODULES([dune-localfunctions], [finiteelements/common/localbasis.hh])
])
