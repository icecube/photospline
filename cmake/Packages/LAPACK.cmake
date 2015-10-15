################################################################################
# Module to find cfitsio                                                       #
#                                                                              #
# This module defines:                                                         #
#                                                                              #
#   LAPACK_FOUND                                                               #
#   LAPACK_VERSION                                                             #
#   LAPACK_LIBRARIES                                                           #
#   LAPACK_INCLUDE_DIR                                                         #
#   LAPACK_LIB_DIR                                                             #
#   LAPACK_CPPFLAGS                                                            #
#   LAPACK_LDFLAGS                                                             #
################################################################################

LIST (APPEND _blas_vendors "Goto" "ATLAS" "Apple" "Generic")

IF (NOT LAPACK_FOUND)

  FOREACH (_vendor ${_blas_vendors})
    SET (BLA_VENDOR ${_vendor})
    FIND_PACKAGE (LAPACK QUIET)
    IF (LAPACK_FOUND)
      SET (LAPACK_LDFLAGS "${LAPACK_LIBRARIES}")
      MESSAGE (STATUS "LAPACK found (${_vendor})")
      MESSAGE (STATUS "  * libs:     ${LAPACK_LIBRARIES}")
      BREAK ()
    ENDIF ()
  ENDFOREACH ()

ENDIF (NOT LAPACK_FOUND)

