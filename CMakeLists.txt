cmake_minimum_required (VERSION 3.1.0 FATAL_ERROR)
cmake_policy(VERSION 3.1.0)

project (photospline VERSION 2.3.0 LANGUAGES C CXX)

SET(CMAKE_CXX_STANDARD 11)
SET(CMAKE_C_STANDARD 99)

# -----------------------------------------------------------------------------
# Set up build paths and external packages
# -----------------------------------------------------------------------------
SET (CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake
                       ${PROJECT_SOURCE_DIR}/cmake/Packages)

INCLUDE (CFITSIO)
INCLUDE (LAPACK)
INCLUDE (SuiteSparse)
INCLUDE (Python)
INCLUDE (package)

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
  ${PROJECT_SOURCE_DIR}/src/core/bspline.cpp
  ${PROJECT_SOURCE_DIR}/src/core/bspline_multi.cpp
  ${PROJECT_SOURCE_DIR}/src/core/convolve.cpp
  ${PROJECT_SOURCE_DIR}/src/core/fitsio.cpp
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
target_compile_options (photospline PRIVATE -Wall -Wextra)
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
install (DIRECTORY ${PROJECT_SOURCE_DIR}/include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# C interface
LIST (APPEND c_HEADERS
  ${PROJECT_SOURCE_DIR}/include/photospline/cinter/splinetable.h
)
LIST (APPEND c_SOURCES
  ${PROJECT_SOURCE_DIR}/src/cinter/splinetable.cpp
)
add_library (cphotospline SHARED
  ${PROJECT_SOURCE_DIR}/src/cinter/splinetable.cpp
)
target_include_directories (cphotospline
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)
target_compile_options (cphotospline PRIVATE -Wall -Wextra)

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
    ${PROJECT_SOURCE_DIR}/src/fitter/cholesky_solve.c
    ${PROJECT_SOURCE_DIR}/src/fitter/glam.c
    ${PROJECT_SOURCE_DIR}/src/fitter/nnls.c
    ${PROJECT_SOURCE_DIR}/src/fitter/splineutil.c
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
    PRIVATE
      ${BLAS_LIBRARIES}
      ${LAPACK_LIBRARIES}
  )

  if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    target_link_libraries(spglam PUBLIC pthread rt)
  endif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  # Assume that if TBB can be found, we need to link
  # against it to use SuiteSparse
  IF (TBB_FOUND)
    target_link_libraries(spglam PRIVATE ${TBB_LIBRARIES})
  ENDIF (TBB_FOUND)
  # Likewise assume that if METIS is found we need it as well, 
  # even if we never use it. 
  IF (METIS_FOUND)
    target_link_libraries(spglam PRIVATE ${METIS_LIBRARIES})
  ENDIF (METIS_FOUND)
  INSTALL (TARGETS spglam EXPORT ${PROJECT_NAME}Config LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
  
  target_compile_definitions(photospline PUBLIC PHOTOSPLINE_INCLUDES_SPGLAM)
  target_link_libraries(photospline PUBLIC spglam)
  target_include_directories(photospline INTERFACE ${SUITESPARSE_INCLUDE_DIR})

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
    ${PROJECT_SOURCE_DIR}/src/python/photosplinemodule.cpp
  )

  add_library (pyphotospline SHARED ${python_SOURCES})
  set_target_properties(pyphotospline PROPERTIES LINKER_LANGUAGE C OUTPUT_NAME "photospline" PREFIX "" SUFFIX ".so")
  target_include_directories (pyphotospline PRIVATE ${PYTHON_INCLUDE_DIRS})
  target_link_libraries (pyphotospline PUBLIC PRIVATE photospline)
  IF(NUMPY_FOUND)
    target_include_directories (pyphotospline PRIVATE ${NUMPY_INCLUDE_DIR})
    target_compile_definitions (pyphotospline PRIVATE HAVE_NUMPY)
    target_compile_definitions (pyphotospline PRIVATE NUMPY_VERSION_MAJOR=${NUMPY_VERSION_MAJOR})
    target_compile_definitions (pyphotospline PRIVATE NUMPY_VERSION_MINOR=${NUMPY_VERSION_MINOR})
    target_compile_definitions (pyphotospline PRIVATE NUMPY_VERSION_PATCH=${NUMPY_VERSION_PATCH})
  ENDIF(NUMPY_FOUND)
  # Assume that undef symbols (e.g. from libpython) will be present at runtime
  if(APPLE)
    SET_TARGET_PROPERTIES(pyphotospline PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
  endif(APPLE)
  EXECUTE_PROCESS (COMMAND ${PYTHON_EXECUTABLE} -c
"from sys import version_info;
if version_info.major < 3 or (version_info.major == 3 and version_info.minor < 10):
 from distutils.sysconfig import get_python_lib; print(get_python_lib(prefix='')); exit(0);
from sysconfig import get_path, get_preferred_scheme;
from pathlib import PurePath;
ps = get_preferred_scheme('prefix');
d  = PurePath(get_path('data', ps));
pl = PurePath(get_path('platlib', ps));
print(pl.relative_to(d))"
    OUTPUT_VARIABLE PYTHON_MODULE_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE)
  INSTALL(TARGETS pyphotospline LIBRARY DESTINATION ${PYTHON_MODULE_DIR})
  IF(${Python_VERSION} VERSION_GREATER_EQUAL "3.9")
    MESSAGE("-- + typing stubs")
    INSTALL(FILES typings/photospline-stubs/__init__.pyi DESTINATION ${PYTHON_MODULE_DIR}/photospline-stubs)
  ENDIF()
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
SET(CONFIG_CFLAGS "${CONFIG_COMPILE_DEFINITIONS} -I${CFITSIO_INCLUDE_DIR}")
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
  SET(CONFIG_CFLAGS "${CONFIG_CFLAGS} -I${SUITESPARSE_INCLUDE_DIR}")
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
CONFIGURE_FILE (${PROJECT_SOURCE_DIR}/cmake/photospline-config.in
  ${PROJECT_BINARY_DIR}/photospline-config
  @ONLY)
INSTALL (PROGRAMS ${PROJECT_BINARY_DIR}/photospline-config DESTINATION bin)

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
  "${PROJECT_SOURCE_DIR}/cmake/cmake_uninstall.cmake.in"
  "${PROJECT_BINARY_DIR}/cmake/cmake_uninstall_cmake"
  IMMEDIATE @ONLY)
ADD_CUSTOM_TARGET (uninstall
  "${CMAKE_COMMAND}" -P "${PROJECT_BINARY_DIR}/cmake/cmake_uninstall_cmake")

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
LIST(APPEND TEST_SOURCES
  test/test_main.cpp
  test/test_fitsio.cpp
  test/test_eval.cpp
)
if(BUILD_SPGLAM)
  LIST(APPEND TEST_SOURCES
    test/test_operators.cpp
  )
endif(BUILD_SPGLAM)

ADD_EXECUTABLE(photospline-test
  ${TEST_SOURCES}
)
TARGET_LINK_LIBRARIES(photospline-test
  photospline
)
target_compile_definitions(photospline-test PRIVATE PHOTOSPLINE_NO_EVAL_TEMPLATES)
target_compile_options (photospline-test PRIVATE -Wsign-compare)
ADD_TEST(photospline-test photospline-test
 WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/test
)
ADD_EXECUTABLE(photospline-test-templated
  ${TEST_SOURCES}
)
TARGET_LINK_LIBRARIES(photospline-test-templated
  photospline
)
target_compile_options (photospline-test-templated PRIVATE -Wsign-compare)
ADD_TEST(photospline-test-templated photospline-test-templated
 WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/test
)
LIST (APPEND ALL_TESTS
  photospline-test
  photospline-test-templated
)

if(PYTHON_FOUND AND NUMPY_FOUND)
  ADD_TEST(photospline-test-pystack ${PYTHON_EXECUTABLE} ${PROJECT_SOURCE_DIR}/test/test_stack.py)
  set_property(TEST photospline-test-pystack PROPERTY ENVIRONMENT PYTHONPATH=${PROJECT_BINARY_DIR})
  LIST (APPEND ALL_TESTS photospline-test-pystack)
  ADD_TEST(photospline-test-pyeval ${PYTHON_EXECUTABLE} ${PROJECT_SOURCE_DIR}/test/test_eval.py)
  set_property(TEST photospline-test-pyeval PROPERTY ENVIRONMENT PYTHONPATH=${PROJECT_BINARY_DIR})
  LIST (APPEND ALL_TESTS photospline-test-pyeval)
endif()

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
   WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/test
  )
  LIST (APPEND ALL_TESTS photospline-test-fit)
  if(NUMPY_FOUND)
    ADD_TEST(photospline-test-pyfit ${PYTHON_EXECUTABLE} ${PROJECT_SOURCE_DIR}/test/test_fit.py
     WORKING_DIRECTORY ${CMAKE_BUILD_DIR}
    )
    LIST (APPEND ALL_TESTS photospline-test-pyfit)
    if(BUILD_SPGLAM)
      ADD_TEST(NAME photospline-test-nnls
        COMMAND ${PYTHON_EXECUTABLE} ${PROJECT_SOURCE_DIR}/test/test_nnls.py
        WORKING_DIRECTORY ${CMAKE_BUILD_DIR}
      )
      set_property(TEST photospline-test-nnls PROPERTY ENVIRONMENT PYTHONPATH=${PROJECT_BINARY_DIR})
      LIST (APPEND ALL_TESTS photospline-test-nnls)
    endif()
  endif()
endif()

ADD_CUSTOM_TARGET(check 
  COMMAND ${CMAKE_CTEST_COMMAND}
  DEPENDS ${ALL_TESTS})

# -----------------------------------------------------------------------------
# Documentation
# -----------------------------------------------------------------------------
ADD_SUBDIRECTORY(docs)
