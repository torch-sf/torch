AC_CONFIG_MACRO_DIRS([m4])

AC_INIT([torch-flash], [0.0])

# Find the compilers
AC_PROG_FC()
AX_FC_ISO_C_BINDING()
AC_PROG_CC()
AC_PROG_CXX()

# Find external dependencies
AC_LANG_PUSH([C])
AX_MPI()
AC_OPENMP()
AC_SUBST(OPENMP_CFLAGS)
AC_LANG_POP([C])


AC_LANG_PUSH([C++])
AX_MPI()
AC_LANG_POP([C++])


AC_LANG_PUSH([Fortran])
AX_MPI()
AC_OPENMP()
AC_SUBST(OPENMP_FCFLAGS)
TORCH_PETSC(3, 19, 0, 3, 20, 0)
AC_CHECK_LIB([HYPRE], [HYPRE_SStructGridDestroy], [LIB_HYPRE=-lHYPRE], [])
AC_SUBST([LIB_HYPRE])
AC_LANG_POP([Fortran])

# This searches all languages
AMUSE_LIB_HDF5()

# Generate output
AC_CONFIG_FILES([Makefile.h])
AC_OUTPUT

