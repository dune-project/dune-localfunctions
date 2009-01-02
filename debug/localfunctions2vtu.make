# -*- mode: makefile -*-

include Makefile

$(PROG_NAME).cc: make-functions2vtu.pl
	perl make-functions2vtu.pl --source "$(BASIS)"

$(PROG_NAME): $(PROG_NAME).o
	@rm -f $(PROG_NAME)
	$(CXXLINK) $(PROG_NAME).o $(LDADD) $(LIBS)

include ./$(DEPDIR)/$(PROG_NAME).Po
