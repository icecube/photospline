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

LIST (APPEND _blas_vendors "OpenBLAS" "Goto" "ATLAS" "Apple" "Generic")

# old CMake doesn't know that OpenBLAS is a BLAS
IF (${CMAKE_VERSION} VERSION_LESS 3.6)
  IF (NOT LAPACK_FOUND)
    set(_vendor OpenBLAS)
    FIND_PACKAGE(${_vendor} QUIET)
    IF (OpenBLAS_FOUND)
      SET (LAPACK_LIBRARIES "${${_vendor}_LIBRARIES}")
      MESSAGE (STATUS "LAPACK found (${_vendor})")
      MESSAGE (STATUS "  * libs:     ${${_vendor}_LIBRARIES}")
      SET (LAPACK_FOUND TRUE)
      SET (BLAS_FOUND TRUE)
    ENDIF ()
  ENDIF ()
ENDIF ()

IF (NOT LAPACK_FOUND)

  FOREACH (_vendor ${_blas_vendors})
    SET (BLA_VENDOR ${_vendor})
    FIND_PACKAGE (LAPACK QUIET)
    IF (LAPACK_FOUND)
      MESSAGE (STATUS "LAPACK found (${_vendor})")
      MESSAGE (STATUS "  * libs:     ${LAPACK_LIBRARIES}")
      BREAK ()
    ENDIF ()
  ENDFOREACH ()

ENDIF (NOT LAPACK_FOUND)

