AC_DEFUN([DUNE_VIRTUAL_BASIS],[
  # enable dynamic polymorphism for shapefunctions
  AC_ARG_ENABLE(virtualshapefunctions,
   AC_HELP_STRING([--enable-virtualshapefunctions],[enable dynamic polymorphism for shapefunctions]))
  if test x$enable_virtualshapefunctions = xyes; then
    AC_DEFINE([DUNE_VIRTUAL_SHAPEFUNCTIONS], [1], 
      [Define to 1 if the shapefunctions should have a virtual basis.])
  fi
  AM_CONDITIONAL([VIRTUAL_SHAPEFUNCTIONS],
    [test x$enable_virtualshapefunctions = xyes])
])