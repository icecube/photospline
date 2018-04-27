cmake_minimum_required (VERSION 3.1.0 FATAL_ERROR)
cmake_policy(VERSION 3.1.0)

project (photospline VERSION 2.0.1 LANGUAGES C CXX)

SET(CMAKE_CXX_STANDARD 11)
SET(CMAKE_C_STANDARD 99)

# -----------------------------------------------------------------------------
# Set up build paths and external packages
# -----------------------------------------------------------------------------
SET (CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake
                       ${CMAKE_SOURCE_DIR}/cmake/Packages)

INCLUDE (CFITSIO)
INCLUDE (LAPACK)
INCLUDE (SuiteSparse)
INCLUDE (Python)

IF (BLAS_FOUND AND LAPACK_FOUND AND SUITESPARSE_FOUND AND NOT DEFINED BUILD_SPGLAM)
  SET(BUILD_SPGLAM TRUE)
ENDIF ()

# -----------------------------------------------------------------------------
# Set up the installation
# -----------------------------------------------------------------------------
# Install to correct subdirectories on all platforms
include(GNUInstallDirs)
if(APPLE)
  # See: https://cmake.org/Wiki/CMake_RPATH_handling
  set(CMAKE_MACOSX_RPATH TRUE)
  set(CMAKE_SKIP_BUILD_RPATH FALSE)
  set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
  set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
endif(APPLE)

# Core library
LIST (APPEND core_SOURCES
  ${CMAKE_SOURCE_DIR}/src/core/bspline.cpp
  ${CMAKE_SOURCE_DIR}/src/core/bspline_multi.cpp
  ${CMAKE_SOURCE_DIR}/src/core/convolve.cpp
  ${CMAKE_SOURCE_DIR}/src/core/fitsio.cpp
)
add_library (photospline SHARED ${core_SOURCES})
target_include_directories (photospline
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
  INTERFACE
    ${CFITSIO_INCLUDE_DIR}
)
target_compile_features (photospline
  PUBLIC
    cxx_constexpr
)
target_compile_options (photospline PUBLIC -O3)
IF (CMAKE_SYSTEM_PROCESSOR MATCHES "^(i.86|x86_64)$")
  target_compile_options (photospline PUBLIC -msse2 -msse3 -msse4 -msse4.1 -msse4.2 -mno-avx)
ELSEIF (CMAKE_SYSTEM_PROCESSOR MATCHES "^(powerpc|ppc)")
  target_compile_options (photospline PUBLIC -maltivec)
ELSEIF (CMAKE_SYSTEM_PROCESSOR MATCHES "^sparc")
  target_compile_options (photospline PUBLIC -mvis)
ENDIF ()
target_link_libraries (photospline
  PUBLIC
    ${CFITSIO_LIBRARIES}
)
target_compile_definitions (photospline
  PUBLIC
    PHOTOSPLINE_VERSION=${PROJECT_VERSION}
)
set_target_properties (photospline
  PROPERTIES 
    SOVERSION ${PROJECT_VERSION_MAJOR}
    VERSION "${PROJECT_VERSION}"
)
install (TARGETS photospline EXPORT ${PROJECT_NAME}Config LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
install (DIRECTORY ${CMAKE_SOURCE_DIR}/include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# C interface
LIST (APPEND c_HEADERS
  ${CMAKE_SOURCE_DIR}/include/photospline/cinter/splinetable.h
)
LIST (APPEND c_SOURCES
  ${CMAKE_SOURCE_DIR}/src/cinter/splinetable.cpp
)
add_library (cphotospline SHARED
  ${CMAKE_SOURCE_DIR}/src/cinter/splinetable.cpp
)
target_include_directories (cphotospline
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

TARGET_LINK_LIBRARIES (cphotospline PUBLIC photospline)
SET_TARGET_PROPERTIES(cphotospline PROPERTIES
  SOVERSION ${PROJECT_VERSION_MAJOR}
  VERSION "${PROJECT_VERSION}"
)
INSTALL (TARGETS cphotospline EXPORT ${PROJECT_NAME}Config LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

# C GLAM fitting library
IF (BUILD_SPGLAM)
  MESSAGE("-- Will build the spglam fitting library")
 
  LIST (APPEND fitter_SOURCES
    ${CMAKE_SOURCE_DIR}/src/fitter/cholesky_solve.c
    ${CMAKE_SOURCE_DIR}/src/fitter/glam.c
    ${CMAKE_SOURCE_DIR}/src/fitter/nnls.c
    ${CMAKE_SOURCE_DIR}/src/fitter/splineutil.c
  )

  ADD_LIBRARY (spglam SHARED ${fitter_SOURCES})
  SET_TARGET_PROPERTIES(spglam PROPERTIES
    SOVERSION ${PROJECT_VERSION_MAJOR}
    VERSION "${PROJECT_VERSION}"
    LINKER_LANGUAGE C
  )
  target_include_directories (spglam
    PUBLIC
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
      $<INSTALL_INTERFACE:include>
    PRIVATE
      ${LAPACK_INCLUDE_DIR}
      ${SUITESPARSE_INCLUDE_DIR}
  )
  target_link_libraries(spglam
    PUBLIC
      ${SUITESPARSE_LIBRARIES}
      ${BLAS_LIBRARIES}
      ${LAPACK_LIBRARIES}
  )

  if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    target_link_libraries(spglam PUBLIC pthread rt)
  endif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  # Assume that if TBB can be found, we need to link
  # against it to use SuiteSparse
  IF (TBB_FOUND)
    target_link_libraries(spglam PUBLIC ${TBB_LIBRARIES})
  ENDIF (TBB_FOUND)
  # Likewise assume that if METIS is found we need it as well, 
  # even if we never use it. 
  IF (METIS_FOUND)
    target_link_libraries(spglam PUBLIC ${METIS_LIBRARIES})
  ENDIF (METIS_FOUND)
  INSTALL (TARGETS spglam EXPORT ${PROJECT_NAME}Config LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
  
  target_compile_definitions(photospline PUBLIC PHOTOSPLINE_INCLUDES_SPGLAM)
  target_link_libraries(photospline PUBLIC spglam)

ELSEIF(NOT BLAS_FOUND)
  MESSAGE("No BLAS library found, skipping build of fitter")
ELSEIF(NOT LAPACK_FOUND)
  MESSAGE("No LAPACK library found, skipping build of fitter")
ELSEIF(NOT SUITESPARSE_FOUND)
  MESSAGE("SuiteSparse not found, skipping build of fitter")
ENDIF()

# Python interface

IF(PYTHON_FOUND)
  MESSAGE("-- Will build the python module")

  LIST (APPEND python_SOURCES
    ${CMAKE_SOURCE_DIR}/src/python/photosplinemodule.cpp
  )

  add_library (pyphotospline SHARED ${python_SOURCES})
  set_target_properties(pyphotospline PROPERTIES LINKER_LANGUAGE C OUTPUT_NAME "photospline" PREFIX "" SUFFIX ".so")
  target_include_directories (pyphotospline PRIVATE ${PYTHON_INCLUDE_DIR})
  target_link_libraries (pyphotospline PRIVATE photospline)
  IF(NUMPY_FOUND)
    target_include_directories (pyphotospline PRIVATE ${NUMPY_INCLUDE_DIR})
    target_compile_definitions (pyphotospline PRIVATE HAVE_NUMPY)
  ENDIF(NUMPY_FOUND)
  # Assume that undef symbols (e.g. from libpython) will be present at runtime
  if(APPLE)
    SET_TARGET_PROPERTIES(pyphotospline PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
  endif(APPLE)
  EXECUTE_PROCESS (COMMAND ${PYTHON_EXECUTABLE} -c
    "from distutils.sysconfig import get_python_lib; print get_python_lib(prefix='')"
    OUTPUT_VARIABLE PYTHON_MODULE_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE)
  INSTALL(TARGETS pyphotospline LIBRARY DESTINATION ${PYTHON_MODULE_DIR})
ELSEIF(NOT PYTHON_FOUND)
  MESSAGE("-- Python not found, skipping build of python module")
ENDIF(PYTHON_FOUND)

# -----------------------------------------------------------------------------
# Create a script for configuring the user environment
# -----------------------------------------------------------------------------
MESSAGE ("-- Creating photospline-config script")
#these may not actually be suitable for users
#get_target_property(CONFIG_COMPILE_OPTIONS_RAW photospline COMPILE_OPTIONS)
#foreach(flag IN LISTS CONFIG_COMPILE_OPTIONS_RAW)
#  SET(CONFIG_COMPILE_OPTIONS "${CONFIG_COMPILE_OPTIONS} ${flag}")
#endforeach(flag)
GET_TARGET_PROPERTY(CONFIG_COMPILE_DEFINITIONS_RAW photospline COMPILE_DEFINITIONS)
FOREACH(def IN LISTS CONFIG_COMPILE_DEFINITIONS_RAW)
  SET(CONFIG_COMPILE_DEFINITIONS "${CONFIG_COMPILE_DEFINITIONS} -D${def}")
ENDFOREACH(def)
SET(CONFIG_CFLAGS "${CONFIG_COMPILE_DEFINITIONS}")
SET(CONFIG_LDFLAGS "-lphotospline")
GET_TARGET_PROPERTY(CONFIG_LINK_LIBRARIES_RAW photospline LINK_LIBRARIES)
FOREACH(lib IN LISTS CONFIG_LINK_LIBRARIES_RAW)
  IF("${lib}" MATCHES ".*lib.*\.(a|so|dylib)$")
    SET(CONFIG_LDFLAGS "${CONFIG_LDFLAGS} ${lib}")
  ELSEIF("${lib}" MATCHES ".*\.framework")
    #ignore
  ELSE()
    SET(CONFIG_LDFLAGS "${CONFIG_LDFLAGS} -l${lib}")
  ENDIF()
ENDFOREACH(lib)
IF(BUILD_SPGLAM)
  GET_TARGET_PROPERTY(CONFIG_LINK_LIBRARIES_RAW spglam LINK_LIBRARIES)
  FOREACH(lib IN LISTS CONFIG_LINK_LIBRARIES_RAW)
    IF("${lib}" MATCHES ".*lib.*\.(a|so|dylib)$")
	  SET(CONFIG_LDFLAGS "${CONFIG_LDFLAGS} ${lib}")
	ELSEIF("${lib}" MATCHES ".*\.framework")
	  #ignore
	ELSE()
	  SET(CONFIG_LDFLAGS "${CONFIG_LDFLAGS} -l${lib}")
	ENDIF()
  ENDFOREACH(lib)
ENDIF(BUILD_SPGLAM)
CONFIGURE_FILE (${CMAKE_SOURCE_DIR}/cmake/photospline-config.in
  ${CMAKE_BINARY_DIR}/photospline-config
  @ONLY)
INSTALL (PROGRAMS ${CMAKE_BINARY_DIR}/photospline-config DESTINATION bin)

# -----------------------------------------------------------------------------
# Export targets for use in downstream CMake projects
# -----------------------------------------------------------------------------
include(CMakePackageConfigHelpers)
# Make importable from build directory
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY AnyNewerVersion
)
export(EXPORT ${PROJECT_NAME}Config FILE ${PROJECT_NAME}Config.cmake)

# Make importable from install location
set(_config_dir share/${PROJECT_NAME}/cmake)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
  DESTINATION ${_config_dir}
)
install(EXPORT ${PROJECT_NAME}Config
  DESTINATION ${_config_dir}
)

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
TARGET_LINK_LIBRARIES(photospline-inspect
  photospline
)
install(TARGETS photospline-inspect RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

ADD_EXECUTABLE(photospline-bench
  src/tools/bench.cpp
)
target_compile_definitions(photospline-bench PRIVATE PHOTOSPLINE_NO_EVAL_TEMPLATES)
TARGET_LINK_LIBRARIES(photospline-bench
  photospline
)
ADD_EXECUTABLE(photospline-bench-templated
  src/tools/bench.cpp
)
TARGET_LINK_LIBRARIES(photospline-bench-templated
  photospline
)

ADD_EXECUTABLE(photospline-eval
  src/tools/eval.cpp
)
TARGET_LINK_LIBRARIES(photospline-eval
  photospline
)
install(TARGETS photospline-eval RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

if(BUILD_SPGLAM)
  ADD_EXECUTABLE(photospline-gen_test_splines
    src/tools/gen_test_splines.cpp
  )
  TARGET_LINK_LIBRARIES(photospline-gen_test_splines
    photospline
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
TARGET_LINK_LIBRARIES(photospline-test
  photospline
)
target_compile_definitions(photospline-test PRIVATE PHOTOSPLINE_NO_EVAL_TEMPLATES)
ADD_TEST(photospline-test photospline-test
 WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/test
)
ADD_EXECUTABLE(photospline-test-templated
  test/test_main.cpp
  test/test_fitsio.cpp
  test/test_eval.cpp
)
TARGET_LINK_LIBRARIES(photospline-test-templated
  photospline
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
  TARGET_LINK_LIBRARIES(photospline-test-fit
    photospline
    spglam
  )
  ADD_TEST(photospline-test-fit photospline-test-fit
   WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/test
  )
  LIST (APPEND ALL_TESTS photospline-test-fit)
  if(NUMPY_FOUND)
    ADD_TEST(photospline-test-pyfit ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/test/test_fit.py
     WORKING_DIRECTORY ${CMAKE_BUILD_DIR}
    )
    LIST (APPEND ALL_TESTS photospline-test-pyfit)
  endif()
endif()

ADD_CUSTOM_TARGET(check 
  COMMAND ${CMAKE_CTEST_COMMAND}
  DEPENDS ${ALL_TESTS})

# -----------------------------------------------------------------------------
# Documentation
# -----------------------------------------------------------------------------
ADD_SUBDIRECTORY(docs)
