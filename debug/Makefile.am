dist_EXTRA_SCRIPTS = make-functions2vtu.pl localfunctions2vtu.make

#use like this: ./make-localfunctions2vtu.pl 'P0LocalBasis<double, double, 2>'
#make-localfunctions2vtu.pl will call the makefile with the two variables
#PROG_NAME and BASIS set on the command line and with the apropriate target

noinst_HEADERS = localfunctions2vtu.hh

EXTRA_PROGRAMS = dummy

dummy_SOURCES = dummy.cc

CLEANFILES = functions2vtu-*

include $(top_srcdir)/am/global-rules

