localfunctionsdir = $(includedir)/dune/localfunctions

localfunctions_HEADERS = lagrange.hh orthonormal.hh raviartthomas.hh \
        edger02d.hh edges02d.hh edges03d.hh \
        p0.hh p1.hh pk2d.hh pk3d.hh q22d.hh prismp1.hh prismp2.hh \
	rannacher_turek2d.hh hierarchicalprismp2.hh \
        monom.hh p2.hh mimetic.hh \
	pk2d.hh q22d.hh q1.hh \
	pq22d.hh

SUBDIRS = common utility test \
          lagrange orthonormal \
          p0 p1 pk2d pk3d q1 q22d raviartthomas edger02d edges02d \
 	  edges03d rannacher_turek2d monom mimetic \
          prismp1 prismp2 hierarchicalprismp2 refined

include $(top_srcdir)/am/global-rules
