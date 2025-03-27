
AC_PREREQ([2.69])
AC_CONFIG_MACRO_DIRS([m4 /usr/share/aclocal])
AC_INIT([Primes], [1.0], [https://github.com/JekasObps])
AC_CONFIG_AUX_DIR([build-aux])


AX_VALGRIND_DFLT([memcheck], [on])
AX_VALGRIND_DFLT([helgrind], [off])
AX_VALGRIND_DFLT([sgcheck], [off])
AX_VALGRIND_DFLT([drd], [off])
AX_VALGRIND_CHECK

AC_CONFIG_SRCDIR([primes.c])

AM_INIT_AUTOMAKE([-Wall -Werror foreign 1.11.2])

AM_PROG_AR
LT_INIT
AC_PROG_CC

# Checks for libraries.
PKG_CHECK_MODULES([CHECK], [check >= 0.9.6])

# Checks for header files.
AC_CHECK_HEADERS([stddef.h stdlib.h string.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_TYPE_SIZE_T
AC_TYPE_SSIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([memmove memset])

# Output files
AC_CONFIG_HEADERS([config.h])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
