SUBDIRS = doxygen

TEXSOURCES = dune-localfunctions-manual.tex

if BUILD_DOCS
  DOCFILES = dune-localfunctions-manual.pdf
  EXTRA_DIST = CMakeLists.txt $(DOCFILES)
  EXTRAINSTALL = $(DOCFILES)
endif

# setting like in dune-web
BASEDIR=../..
CURDIR=doc/localfunctions
docdir=$(datadir)/doc/dune-localfunctions

include $(top_srcdir)/am/global-rules
include $(top_srcdir)/am/webstuff
include $(top_srcdir)/am/latex
