INCLUDES = @GSL_CFLAGS@ -I $(top_srcdir)/include/
lib_LTLIBRARIES = libcspl.la
libcspl_la_includedir = $(top_srcdir)/include/cspl
libcspl_la_include_HEADRS = cspl.h
libcspl_la_SOURCES = cspl_math.c cspl_qrs.c cspl_qrs_fit.c 
libcspl_la_LIBADD = @GSL_LIBS@
pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = cspl.pc

