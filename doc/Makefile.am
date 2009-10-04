# $Id: Makefile.am 836 2007-12-05 17:34:43Z robertk $

SUBDIRS = doxygen

EXTRA_DIST = Makefile.dist.in Makefile.dist.am

# setting like in dune-web
BASEDIR=../..
CURDIR=doc/localfunctions
docdir=$(datadir)/doc/dune-localfunctions

doc_DATA = $(DOCS)

# files to clean only in svn tree but not in tarball
SVNCLEANFILES =  *.pdf
clean-local:
	if test -e $(top_srcdir)/doc/doxygen/Doxydep; then rm -rf $(SVNCLEANFILES); fi

if BUILD_DOCS
dist-hook:
	sed $(srcdir)/Makefile.dist.am -e 's/Makefile\.dist/Makefile/g' > $(distdir)/Makefile.am
	sed $(srcdir)/Makefile.dist.in -e 's/Makefile\.dist/Makefile/g' > $(distdir)/Makefile.in
else
dist-hook:
	echo "# No documentation included in distribution! " > $(distdir)/$(DOCUMENTATION_TAG_FILE)
	sed $(srcdir)/Makefile.dist.am -e 's/Makefile\.dist/Makefile/g' > $(distdir)/Makefile.am
	sed $(srcdir)/Makefile.dist.in -e 's/Makefile\.dist/Makefile/g' > $(distdir)/Makefile.in
endif

include $(top_srcdir)/am/latex
include $(top_srcdir)/am/webstuff
include $(top_srcdir)/am/global-rules
