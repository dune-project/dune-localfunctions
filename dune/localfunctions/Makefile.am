localfunctionsdir = $(includedir)/dune/localfunctions

localfunctions_HEADERS = lagrange.hh orthonormal.hh raviartthomas.hh rannacherturek.hh \
        monom.hh mimetic.hh \
	refined.hh

SUBDIRS = common utility test \
          lagrange orthonormal \
	mimetic					\
	monom					\
	rannacherturek			\
        raviartthomas				\
          whitney                               \
          hierarchical refined

include $(top_srcdir)/am/global-rules
