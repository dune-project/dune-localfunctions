localfunctionsdir = $(includedir)/dune/localfunctions

localfunctions_HEADERS = lagrange.hh orthonormal.hh raviartthomas.hh \
        edges02d.hh					\
	edges03d.hh					\
	hierarchicalprismp2.hh \
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
          hierarchicalprismp2 refined

include $(top_srcdir)/am/global-rules
