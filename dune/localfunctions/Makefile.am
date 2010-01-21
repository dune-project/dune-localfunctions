localfunctionsdir = $(includedir)/dune/localfunctions

localfunctions_HEADERS = lagrange.hh orthonormal.hh raviartthomas.hh rannacherturek.hh \
        edges02d.hh					\
	edges03d.hh					\
        monom.hh mimetic.hh \
	refined.hh

SUBDIRS = common utility test \
          lagrange orthonormal \
	edges02d				\
	edges03d				\
	mimetic					\
	monom					\
	rannacherturek			\
        raviartthomas				\
          hierarchical refined

include $(top_srcdir)/am/global-rules
