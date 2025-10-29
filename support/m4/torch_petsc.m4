dnl PETSc detection for Torch

dnl Detect PETSc availability and whether it is sufficiently new

dnl Usage:
dnl TORCH_PETSC(ge_major, ge_minor, ge_subminor, lt_major, lt_minor, lt_subminor)

dnl where ge_major.ge_minor.ge_subminor is the minimum version needed, and
dnl lt_major.lt_minor.lt_subminor is the minimum version that is no longer accepted. So
dnl TORCH_PETSC(3.20.0, 3.21.0) accepts any 3.20.x version, but not 3.21.0.

dnl Output:
dnl Sets PETSC_FOUND to "yes" if a suitable version of PETSc was found, else to "no".
dnl Sets PETSC_CFLAGS, PETSC_LDFLAGS, PETSC_LIBS, and calls AC_SUBST on them.


AC_DEFUN([TORCH_PETSC], [
    dnl Check for the PETSC_DIR and optionally PETSC_ARCH environment variables that
    dnl should exist if PETSc was installed from source, and will point us to it.

    torch_petsc_ge_major="$1"
    torch_petsc_ge_minor="$2"
    torch_petsc_ge_subminor="$3"
    torch_petsc_ge_version="$1.$2.$3"

    torch_petsc_lt_major="$4"
    torch_petsc_lt_minor="$5"
    torch_petsc_lt_subminor="$6"
    torch_petsc_lt_version="$4.$5.$6"

    AC_LANG_PUSH([C])

    dnl Find PETSc based on PETSC_DIR or on the system
    AC_MSG_CHECKING([for PETSc])

    AS_IF([test -n "${PETSC_DIR}"], [
        dnl Source build, in-tree and out-of-tree
        AS_IF([test -f "${PETSC_DIR}/include/petscversion.h"], [
            torch_petsc_found="yes"
            torch_petsc_location=", at ${PETSC_DIR}/${PETSC_ARCH}"
            torch_petsc_cflags="-I${PETSC_DIR}/include -I${PETSC_DIR}/${PETSC_ARCH}/include"
            torch_petsc_ldflags="-L${PETSC_DIR}/${PETSC_ARCH}/lib"
        ], [
            dnl apt on Debian/Ubuntu with PETSC_DIR set to /usr
            AS_IF([test -f "${PETSC_DIR}/include/petsc/petscversion.h"], [
                torch_petsc_found="yes"
                torch_petsc_location=", at ${PETSC_DIR}"
                torch_petsc_cflags="-I${PETSC_DIR}/include/petsc"
                torch_petsc_ldflags="-L${PETSC_DIR}/lib"
            ])
        ])
    ], [
        dnl With Conda, no flags are needed
        AC_COMPILE_IFELSE([AC_LANG_SOURCE([[
            #include "petscversion.h"

            int main() {}
        ]])], [
            torch_petsc_found="yes"
            torch_petsc_location=", in the environment"
            torch_petsc_cflags=""
            torch_petsc_ldflags=""
        ], [
            dnl apt on Debian/Ubuntu puts it here
            AS_IF([test -f /usr/include/petsc/petscversion.h], [
                torch_petsc_found="yes"
                torch_petsc_location=", at /usr/include"
                torch_petsc_cflags="-I/usr/include/petsc"
                torch_petsc_ldflags=""
            ], [
                torch_petsc_found="no"
                torch_petsc_location=""
                torch_petsc_cflags=""
                torch_petsc_ldflags=""
            ])
        ])
    ])

    AC_MSG_RESULT([${torch_petsc_found}${torch_petsc_location}])

    dnl Verify that we can also link

    AC_MSG_CHECKING([whether we can link with PETSc])

    torch_petsc_save_cflags="$CFLAGS"
    torch_petsc_save_ldflags="$LDFLAGS"
    torch_petsc_save_libs="$LIBS"

    CFLAGS="$CFLAGS $torch_petsc_cflags"
    LDFLAGS="$LDFLAGS $torch_petsc_ldflags"
    LIBS="$LIBS -lpetsc"

    AC_LINK_IFELSE([AC_LANG_SOURCE([[
        #include "petscversion.h"

        int main() {}
    ]])], [], [
        torch_petsc_found="no"
    ])

    LIBS="$torch_petsc_save_libs"
    LDFLAGS="$torch_petsc_save_ldflags"
    CFLAGS="$torch_petsc_save_cflags"

    AC_MSG_RESULT([$torch_petsc_found])


    AS_IF([test "$torch_petsc_found" = "yes"], [

        dnl Check version lower bound

        torch_petsc_save_CFLAGS="$CFLAGS"
        CFLAGS="$torch_petsc_cflags ${CFLAGS}"

        AC_MSG_CHECKING([for PETSc >=$torch_petsc_ge_version])

        AC_COMPILE_IFELSE([AC_LANG_SOURCE([[
            #include "petscversion.h"

            int main() {

            #if PETSC_VERSION_LT($torch_petsc_ge_major, $torch_petsc_ge_minor, $torch_petsc_ge_subminor)
                ERROR
            #endif
            }
        ]])], [
            AC_MSG_RESULT([yes])
            torch_petsc_too_old=no
        ], [
            AC_MSG_RESULT([no])
            torch_petsc_too_old=yes
        ])


        dnl Check version upper bound

        AC_MSG_CHECKING([for PETSc <$torch_petsc_lt_version])

        AC_COMPILE_IFELSE([AC_LANG_SOURCE([[
            #include "petscversion.h"

            int main() {

            #if PETSC_VERSION_GE($torch_petsc_lt_major, $torch_petsc_lt_minor, $torch_petsc_lt_subminor)
                ERROR
            #endif
            }
        ]])], [
            AC_MSG_RESULT([yes])
            torch_petsc_too_new=no
        ], [
            AC_MSG_RESULT([no])
            torch_petsc_too_new=yes
        ])

        CFLAGS="$torch_petsc_save_CFLAGS"

    ])

    AC_LANG_POP([C])

    AS_IF([test "$torch_petsc_too_old" = "no" && test "$torch_petsc_too_new" = "no"], [
        PETSC_FOUND=yes
        PETSC_CFLAGS="$torch_petsc_cflags"
        PETSC_LDFLAGS="$torch_petsc_ldflags"
        PETSC_LIBS='-lpetsc'
    ])

    AC_SUBST(PETSC_CFLAGS)
    AC_SUBST(PETSC_LDFLAGS)
    AC_SUBST(PETSC_LIBS)
])

