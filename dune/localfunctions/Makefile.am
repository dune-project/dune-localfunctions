localfunctionsdir = $(includedir)/dune/localfunctions

localfunctions_HEADERS = lagrange.hh orthonormal.hh raviartthomas.hh \
        edger02d.hh edges02d.hh edges03d.hh \
	rannacher_turek2d.hh hierarchicalprismp2.hh \
        monom.hh mimetic.hh \
	refined.hh

SUBDIRS = common utility test \
          lagrange orthonormal \
          raviartthomas edger02d edges02d \
 	  edges03d rannacher_turek2d monom mimetic \
          hierarchicalprismp2 refined

include $(top_srcdir)/am/global-rules
