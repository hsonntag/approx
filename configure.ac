#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

# AC_PREREQ([2.63])
AC_INIT([QRS Splineapproximation], [0.1], [hermann.sonntag@tu-ilmenau.de])
AM_INIT_AUTOMAKE(approx, 0.1)
LT_INIT([disable-shared])

######################################################################
# DOXYGEN SUPPORT
######################################################################

DX_HTML_FEATURE(ON)
DX_CHM_FEATURE(OFF)
DX_CHI_FEATURE(OFF)
DX_MAN_FEATURE(ON)
DX_RTF_FEATURE(OFF)
DX_XML_FEATURE(OFF)
DX_PDF_FEATURE(ON)
DX_PS_FEATURE(OFF)

DX_INIT_DOXYGEN([$PACKAGE_NAME],[doxygen.cfg], [doc])

######################################################################

AC_CONFIG_SRCDIR([config.h.in])
AC_CONFIG_MACRO_DIR([m4])
#AC_CONFIG_HEADERS([config.h])
AM_CONFIG_HEADER(config.h)

# Checks for programs.
AC_PROG_CC
AC_PROG_LIBTOOL
AC_PROG_CXX
AC_PROG_AWK
AC_PROG_CPP
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET
##AC_PROG_RANLIB


# Checks for libraries.
PKG_CHECK_MODULES(GSL, gsl)
AC_SUBST(GSL_CFLAGS)
AC_SUBST(GSL_LIBS)

PKG_CHECK_MODULES(SP, cspl)
AC_SUBST(SP_CFLAGS)
AC_SUBST(SP_LIBS)
# AC_CHECK_LIB([gsl], [gsl_cblas_dgemm])

#AC_CHECK_LIB([gslcblas], [cblas_dgemm])

AC_CHECK_LIB([m], [sin])

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h string.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
#AC_FUNC_STRNLEN
AC_CHECK_FUNCS([pow sqrt strtol strnlen])

# Try to enforce ansi C (ISO C99) with warnings generation
TMP_CFLAGS="$CFLAGS"
CFLAGS="$CFLAGS -ansi -std=c99 -pedantic -Wall"
AC_MSG_CHECKING([for ISO C99 warnings generation])
AC_TRY_COMPILE(
,
int 
main(void) 
{
   return 0L;
},
AC_MSG_RESULT([yes]),
AC_MSG_RESULT([no])
CFLAGS="$TMP_CFLAGS"
)

AC_CONFIG_FILES([Makefile \
		 lib/Makefile \
		 lib/cspl.pc \
		 qrs/Makefile ])
AC_OUTPUT


