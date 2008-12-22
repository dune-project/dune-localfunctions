noinst_PROGRAMS =

if DUNE_GRID
  noinst_PROGRAMS += dune-localfunctions2vtu
endif DUNE_GRID
dune_localfunctions2vtu_SOURCES = dune-localfunctions2vtu.cc
dune_localfunctions2vtu_CPPFLAGS = $(AM_CPPFLAGS) \
	-DLOCAL_BASIS_HEADER='<dune/finiteelements/p0/p0localbasis.hh>' \
	-DLOCAL_BASIS_TYPE='Dune::P0LocalBasis<double, double, 2>' \
	-DLOCAL_BASIS_TYPE_S='"P0LocalBasis<double, double, 2>"' \
	-DPROG_NAME='"dune-localfunctions2vtu-p0-2d"'

include $(top_srcdir)/am/global-rules

