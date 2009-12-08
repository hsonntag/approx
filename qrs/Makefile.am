INCLUDES = @SP_CFLAGS@ -I$(top_srcdir)/lib/ -I$(top_srcdir)/include/
bin_PROGRAMS = approx_qrs
approx_qrs_SOURCES = main.c $(top_srcdir)/lib/cmdline.c
approx_qrs_LDADD = -L$(top_srcdir)/lib/ @SP_LIBS@
