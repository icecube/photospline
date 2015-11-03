PROJECT (photospline C CXX)

CMAKE_MINIMUM_REQUIRED (VERSION 2.8)

# -----------------------------------------------------------------------------
# Release information.
# -----------------------------------------------------------------------------
SET (PHOTOSPLINE_MAJOR_VERSION 2 CACHE STRING "Major version number")
SET (PHOTOSPLINE_MINOR_VERSION 0 CACHE STRING "Minor version number")
SET (PHOTOSPLINE_PATCH_VERSION 0 CACHE STRING "Patch number")
SET (PHOTOSPLINE_CODENAME "trunk" CACHE STRING "Code/type name for project")

# -----------------------------------------------------------------------------
# Set up build paths and external packages
# -----------------------------------------------------------------------------
SET (CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake
                       ${CMAKE_SOURCE_DIR}/cmake/Packages)

INCLUDE (CFITSIO)
INCLUDE (LAPACK)
INCLUDE (SuiteSparse)
INCLUDE (METIS)
INCLUDE (Python)

IF (BLAS_FOUND AND LAPACK_FOUND AND SUITESPARSE_FOUND AND NOT DEFINED BUILD_SPGLAM)
  SET(BUILD_SPGLAM TRUE)
ENDIF ()

SET (PHOTOSPLINE_BUILD "${CMAKE_BINARY_DIR}")
SET (PHOTOSPLINE_SRC "${CMAKE_SOURCE_DIR}")
SET (PHOTOSPLINE_INSTALL "${CMAKE_INSTALL_PREFIX}")
SET (PHOTOSPLINE_BINDIR "${CMAKE_INSTALL_PREFIX}/bin")
SET (PHOTOSPLINE_LIBDIR "${CMAKE_INSTALL_PREFIX}/lib")
SET (PHOTOSPLINE_INCDIR "${CMAKE_INSTALL_PREFIX}/include")
IF (NOT DEFINED PHOTOSPLINE_CFLAGS)
  SET (PHOTOSPLINE_CFLAGS "-I${CMAKE_SOURCE_DIR}/include -O3 -std=c++11")
  IF (CMAKE_SYSTEM_PROCESSOR MATCHES "^i.86$")
    SET (PHOTOSPLINE_CFLAGS "${PHOTOSPLINE_CFLAGS} -msse2 -mfpmath=sse")
  ELSEIF (CMAKE_SYSTEM_PROCESSOR MATCHES "^(powerpc|ppc)")
    SET (PHOTOSPLINE_CFLAGS "${PHOTOSPLINE_CFLAGS} -maltivec")
  ELSEIF (CMAKE_SYSTEM_PROCESSOR MATCHES "^sparc")
    SET (PHOTOSPLINE_CFLAGS "${PHOTOSPLINE_CFLAGS} -mvis")
  ENDIF ()
  IF (BUILD_SPGLAM)
    SET (PHOTOSPLINE_CFLAGS "${PHOTOSPLINE_CFLAGS} -DPHOTOSPLINE_INCLUDES_SPGLAM")
    STRING(REGEX REPLACE "-std=c\\+\\+11" "" SPGLAM_CFLAGS PHOTOSPLINE_CFLAGS)
    SET (SPGLAM_CFLAGS ${SPGLAM_CFLAGS} CACHE STRING "Compiler flags specific to this project")
  ENDIF ()
  SET (PHOTOSPLINE_CFLAGS ${PHOTOSPLINE_CFLAGS} CACHE STRING "Compiler flags specific to this project")
ENDIF ()
#ADD_DEFINITIONS (${PLATFORM_FLAGS} ${PHOTOSPLINE_CFLAGS})

# -----------------------------------------------------------------------------
# Set up the installation
# -----------------------------------------------------------------------------
INCLUDE_DIRECTORIES (
  ${CMAKE_SOURCE_DIR}/include
  ${CFITSIO_INCLUDE_DIR}
)

# Core library
LIST (APPEND core_HEADERS
  ${CMAKE_SOURCE_DIR}/include/photospline/bspline.h
  ${CMAKE_SOURCE_DIR}/include/photospline/splinetable.h
  ${CMAKE_SOURCE_DIR}/include/photospline/detail/bspline_eval.h
  ${CMAKE_SOURCE_DIR}/include/photospline/detail/bspline_multi.h
  ${CMAKE_SOURCE_DIR}/include/photospline/detail/aux.h
  ${CMAKE_SOURCE_DIR}/include/photospline/detail/convolve.h
  ${CMAKE_SOURCE_DIR}/include/photospline/detail/fit.h
  ${CMAKE_SOURCE_DIR}/include/photospline/detail/fitsio.h
  ${CMAKE_SOURCE_DIR}/include/photospline/detail/simd.h
  ${CMAKE_SOURCE_DIR}/include/photospline/detail/sample.h
)
LIST (APPEND core_SOURCES
  ${CMAKE_SOURCE_DIR}/src/core/bspline.cpp
  ${CMAKE_SOURCE_DIR}/src/core/bspline_multi.cpp
  ${CMAKE_SOURCE_DIR}/src/core/convolve.cpp
  ${CMAKE_SOURCE_DIR}/src/core/fitsio.cpp
)
ADD_LIBRARY (photospline-mod SHARED ${core_SOURCES})
TARGET_LINK_LIBRARIES (photospline-mod
  ${CFITSIO_LIBRARIES}
)
SET_TARGET_PROPERTIES(photospline-mod PROPERTIES 
  COMPILE_FLAGS "${PHOTOSPLINE_CFLAGS}"
)
INSTALL (TARGETS photospline-mod LIBRARY DESTINATION lib)
FOREACH(header ${core_HEADERS})
  FILE(RELATIVE_PATH relpath ${CMAKE_SOURCE_DIR}/include ${header})
  GET_FILENAME_COMPONENT(dir ${relpath} DIRECTORY)
  INSTALL(FILES ${header} DESTINATION ${PHOTOSPLINE_INCDIR}/${dir})
endforeach(header ${core_HEADERS})
SET (PHOTOSPLINE_LDFLAGS "-L${PHOTOSPLINE_LIBDIR} -lphotospline-mod")

SET (PHOTOSPLINE_LDFLAGS ${PHOTOSPLINE_LDFLAGS}
     CACHE STRING "Linker flags specific to this project")

# C interface
LIST (APPEND c_HEADERS
  ${CMAKE_SOURCE_DIR}/include/photospline/cinter/splinetable.h
)
LIST (APPEND c_SOURCES
  ${CMAKE_SOURCE_DIR}/src/cinter/splinetable.cpp
)
ADD_LIBRARY (cphotospline SHARED ${c_SOURCES})
LIST (APPEND cinter_LINK_LIBRARIES
  ${CFITSIO_LIBRARIES}
  photospline-mod
)
IF (BUILD_SPGLAM)
  LIST (APPEND cinter_LINK_LIBRARIES
    spglam
  )
ENDIF ()
TARGET_LINK_LIBRARIES (cphotospline ${cinter_LINK_LIBRARIES})
ADD_DEPENDENCIES(cphotospline photospline-mod)
SET_TARGET_PROPERTIES(cphotospline PROPERTIES 
  COMPILE_FLAGS "${PHOTOSPLINE_CFLAGS} -DPHOTOSPLINE_EVAL_TEMPLATES"
  LINK_FLAGS "-L${CMAKE_BINARY_DIR}"
)
INSTALL (TARGETS cphotospline LIBRARY DESTINATION lib)
FOREACH(header ${c_HEADERS})
  FILE(RELATIVE_PATH relpath ${CMAKE_SOURCE_DIR}/include ${header})
  GET_FILENAME_COMPONENT(dir ${relpath} DIRECTORY)
  INSTALL(FILES ${header} DESTINATION ${PHOTOSPLINE_INCDIR}/${dir})
endforeach(header ${c_HEADERS})

# C GLAM fitting library
IF (BUILD_SPGLAM)
  MESSAGE("-- Will build the spglam fitting library")
  INCLUDE_DIRECTORIES (
    ${LAPACK_INCLUDE_DIR}
    ${SUITESPARSE_INCLUDE_DIR})

  LIST (APPEND fitter_HEADDERS
    ${CMAKE_SOURCE_DIR}/include/detail/glam.h
    ${CMAKE_SOURCE_DIR}/include/detail/splineutil.h
  ) 
 
  LIST (APPEND fitter_SOURCES
    ${CMAKE_SOURCE_DIR}/src/fitter/cholesky_solve.c
    ${CMAKE_SOURCE_DIR}/src/fitter/glam.c
    ${CMAKE_SOURCE_DIR}/src/fitter/nnls.c
    ${CMAKE_SOURCE_DIR}/src/fitter/splineutil.c
  )

  ADD_LIBRARY (spglam SHARED ${fitter_SOURCES})
  LIST (APPEND fitter_LINK_LIBRARIES
    ${CFITSIO_LIBRARIES}
    ${SUITESPARSE_LIBRARIES}
    ${BLAS_LIBRARIES}
    ${LAPACK_LIBRARIES}
    photospline-mod
  )
  SET_TARGET_PROPERTIES(spglam PROPERTIES
    COMPILE_FLAGS "${spglam_CFLAGS}"
  )
  # Assume that if TBB can be found, we need to link
  # against it to use SuiteSparse
  IF (TBB_FOUND)
    LIST (APPEND fitter_LINK_LIBRARIES ${TBB_LIBRARIES})
  ENDIF (TBB_FOUND)
  # Likewsie assume that if METIS is found we need it as well, 
  # even if we never use it. 
  IF (METIS_FOUND)
    LIST (APPEND fitter_LINK_LIBRARIES ${METIS_LIBRARIES})
  ENDIF (METIS_FOUND)
  SET_TARGET_PROPERTIES(spglam PROPERTIES LINKER_LANGUAGE C)
  TARGET_LINK_LIBRARIES (spglam ${fitter_LINK_LIBRARIES})
  INSTALL (TARGETS spglam LIBRARY DESTINATION lib)
  INSTALL (FILES ${fitter_HEADERS} DESTINATION ${PHOTOSPLINE_INCDIR}/photospline)

ELSEIF(NOT BLAS_FOUND)
  MESSAGE("No BLAS library found, skipping build of fitter")
ELSEIF(NOT LAPACK_FOUND)
  MESSAGE("No LAPACK library found, skipping build of fitter")
ELSEIF(NOT SUITESPARSE_FOUND)
  MESSAGE("SuiteSparse not found, skipping build of fitter")
ENDIF()

# -----------------------------------------------------------------------------
# Create a script for configuring the user environment
# -----------------------------------------------------------------------------
MESSAGE ("-- Creating photospline-config script")
CONFIGURE_FILE (${CMAKE_SOURCE_DIR}/cmake/photospline-config.in
  ${CMAKE_BINARY_DIR}/photospline-config
  @ONLY)
INSTALL (PROGRAMS ${CMAKE_BINARY_DIR}/photospline-config DESTINATION bin)

# -----------------------------------------------------------------------------
# Add uninstall target for running "make uninstall" in the build directory
# -----------------------------------------------------------------------------
CONFIGURE_FILE (
  "${CMAKE_SOURCE_DIR}/cmake/cmake_uninstall.cmake.in"
  "${CMAKE_BINARY_DIR}/cmake/cmake_uninstall_cmake"
  IMMEDIATE @ONLY)
ADD_CUSTOM_TARGET (uninstall
  "${CMAKE_COMMAND}" -P "${CMAKE_BINARY_DIR}/cmake/cmake_uninstall_cmake")

# -----------------------------------------------------------------------------
# Tools
# -----------------------------------------------------------------------------
ADD_EXECUTABLE(photospline-inspect
  src/tools/inspect.cpp
)
SET_TARGET_PROPERTIES(photospline-inspect PROPERTIES 
  COMPILE_FLAGS "${PHOTOSPLINE_CFLAGS}"
)
TARGET_LINK_LIBRARIES(photospline-inspect
  photospline-mod
  spglam
  ${CFITSIO_LIBRARIES}
)

if(BUILD_SPGLAM)
  ADD_EXECUTABLE(photospline-gen_test_splines
    src/tools/gen_test_splines.cpp
  )
  SET_TARGET_PROPERTIES(photospline-gen_test_splines PROPERTIES
    COMPILE_FLAGS "${PHOTOSPLINE_CFLAGS}"
  )
  TARGET_LINK_LIBRARIES(photospline-gen_test_splines
    photospline-mod
    spglam
    ${CFITSIO_LIBRARIES}
  )
endif(BUILD_SPGLAM)

# -----------------------------------------------------------------------------
# Tests
# -----------------------------------------------------------------------------
ENABLE_TESTING()
ADD_EXECUTABLE(photospline-test
  test/test_main.cpp
  test/test_fitsio.cpp
  test/test_eval.cpp
)
SET_TARGET_PROPERTIES(photospline-test PROPERTIES 
  COMPILE_FLAGS "${PHOTOSPLINE_CFLAGS}"
)
TARGET_LINK_LIBRARIES(photospline-test
  photospline-mod
  ${CFITSIO_LIBRARIES}
)
ADD_TEST(photospline-test photospline-test
 WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/test
)
ADD_EXECUTABLE(photospline-test-templated
  test/test_main.cpp
  test/test_fitsio.cpp
  test/test_eval.cpp
)
SET_TARGET_PROPERTIES(photospline-test-templated PROPERTIES 
  COMPILE_FLAGS "${PHOTOSPLINE_CFLAGS} -DPHOTOSPLINE_EVAL_TEMPLATES"
)
TARGET_LINK_LIBRARIES(photospline-test-templated
  photospline-mod
  ${CFITSIO_LIBRARIES}
)
ADD_TEST(photospline-test-templated photospline-test-templated
 WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/test
)
LIST (APPEND ALL_TESTS
  photospline-test
  photospline-test-templated
)

if(BUILD_SPGLAM)
  ADD_EXECUTABLE(photospline-test-fit
    test/test_main.cpp
    test/test_fit.cpp
  )
  SET_TARGET_PROPERTIES(photospline-test-fit PROPERTIES 
    COMPILE_FLAGS "${PHOTOSPLINE_CFLAGS}"
  )
  TARGET_LINK_LIBRARIES(photospline-test-fit
    photospline-mod
    spglam
    ${CFITSIO_LIBRARIES}
  )
  ADD_TEST(photospline-test-fit photospline-test-fit
   WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/test
  )
  LIST (APPEND ALL_TESTS photospline-test-fit)
endif()

ADD_CUSTOM_TARGET(check 
  COMMAND ${CMAKE_CTEST_COMMAND}
  DEPENDS ${ALL_TESTS})

# -----------------------------------------------------------------------------
# Documentation
# -----------------------------------------------------------------------------
ADD_SUBDIRECTORY(docs)
