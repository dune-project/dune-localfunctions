#use like this: ./make-localfunctions2vtu.pl 'P0LocalBasis<double, double, 2>'
#make-localfunctions2vtu.pl will call the makefile with the two variables
#PROG_NAME and BASIS set on the command line and with the apropriate target

dist_EXTRA_SCRIPTS = make-functions2vtu.pl
#this header is mentioned nowhere else, mention it here so it gets distributed
noinst_HEADERS = localfunctions2vtu.hh

# dummy target so automake writes the rules and variables for C++ programs
# into the Makefile.in
EXTRA_PROGRAMS = dummy
dummy_SOURCES = dummy.cc

#set to dummy target by default so the inclusion of the dependency file below
#works
PROG_NAME = dummy

#do nothing if PROG_NAME=dummy
$(PROG_NAME).cc: make-functions2vtu.pl
	case "$(PROG_NAME)" in					\
	dummy) touch dummy.cc;;					\
	*) perl make-functions2vtu.pl --source "$(BASIS)";;	\
	esac

#this linking rule needs to be here in full, otherwise make will use the

#"%: %.cc" rule, which is not overwritten by automake, instead of the
# "%: %.o", "%.o: %.cc" chain, which _is_ overwritten by automake
$(PROG_NAME): $(PROG_NAME).o
	@rm -f $(PROG_NAME)
	$(CXXLINK) $(PROG_NAME).o $(LDADD) $(LIBS)

#include the dependencies for this target.  Include the deps for dummy if
#PROG_NAME is not explicitly set.
include ./$(DEPDIR)/$(PROG_NAME).Po

#clean all generated program sources, the compiled programs themselves, and
#their dependencies
CLEANFILES = functions2vtu-* ./$(DEPDIR)/functions2vtu-*

include $(top_srcdir)/am/global-rules

