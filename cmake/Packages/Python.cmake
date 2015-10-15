################################################################################
# Module to find python                                                        #
#                                                                              #
# This module will call the default CMake python modules and define:           #
#                                                                              #
#   PYTHON_FOUND                                                               #
#   PYTHON_LIBRARIES                                                           #
#   PYTHON_INCLUDE_DIR                                                         #
#   PYTHON_CPPFLAGS                                                            #
#   PYTHON_LDFLAGS                                                             #
#   NUMPY_FOUND
#                                                                              #
# plus the other PYTHON_* variables defined by CMake (see CMake documentation) #
################################################################################

IF (NOT PYTHON_FOUND)
  SET (PYTHONROOT $ENV{PYTHONHOME})

  IF (IS_DIRECTORY ${PYTHONROOT})

    FIND_PATH (PYTHON_EXECUTABLE_DIR python
      PATHS ${PYTHONROOT}/bin
      NO_DEFAULT_PATH)

    IF (IS_DIRECTORY ${PYTHON_EXECUTABLE_DIR})
      SET (PYTHON_EXECUTABLE ${PYTHON_EXECUTABLE_DIR}/python)
      
      EXECUTE_PROCESS (COMMAND ${PYTHON_EXECUTABLE} -V
        ERROR_VARIABLE PYTHON_VERSION_STRING
        ERROR_STRIP_TRAILING_WHITESPACE)
      SEPARATE_ARGUMENTS (PYTHON_VERSION_STRING)

      SET (PYTHON_VERSION "")
      SET (PYTHON_MAJOR_VERSION "")
      FOREACH (item ${PYTHON_VERSION_STRING})
        SET (PYTHON_VERSION ${item})
      ENDFOREACH (item ${PYTHON_VERSION_STRING})

      STRING(LENGTH ${PYTHON_VERSION} POS)
      MATH(EXPR POS "${POS} - 2")
      STRING(SUBSTRING ${PYTHON_VERSION} 0 ${POS} PYTHON_MAJOR_VERSION)

      SET (PYTHON_INCLUDE_PATH ${PYTHONROOT}/include/python${PYTHON_MAJOR_VERSION})

      FIND_PATH (PYTHON_INCLUDE_DIR Python.h
        PATHS ${PYTHON_INCLUDE_PATH} 
        NO_DEFAULT_PATH)

      FIND_LIBRARY (PYTHON_LIBRARIES NAMES python2.7
        PATHS ${PYTHONROOT}/lib
        NO_DEFAULT_PATH)

      SET (PYTHON_CPPFLAGS "-I${PYTHON_INCLUDE_DIR}")
      SET (PYTHON_LDFLAGS "${PYTHON_LIBRARIES}")
    ENDIF (IS_DIRECTORY ${PYTHON_EXECUTABLE_DIR})

    MESSAGE (STATUS "Python executable: ${PYTHON_EXECUTABLE}")

    # Search for numpy
    EXECUTE_PROCESS (COMMAND ${PYTHON_EXECUTABLE} -c "import numpy"
      RESULT_VARIABLE NUMPY_FOUND)
    IF (NUMPY_FOUND EQUAL 0)
      SET (NUMPY_FOUND TRUE)
    ELSE ()
      SET (NUMPY_FOUND FALSE)
    ENDIF ()

    IF (NUMPY_FOUND)
      SET (NUMPY_FOUND TRUE CACHE BOOL "Numpy found successfully" FORCE)
      EXECUTE_PROCESS (COMMAND ${PYTHON_EXECUTABLE} -c
        "import numpy; print numpy.get_include()"
        OUTPUT_VARIABLE NUMPY_INCLUDE_DIR
        OUTPUT_STRIP_TRAILING_WHITESPACE)
      SET (NUMPY_INCLUDE_DIR ${NUMPY_INCLUDE_DIR} CACHE STRING "Numpy directory")
    ELSE ()
      SET (NUMPY_FOUND FALSE CACHE BOOL "Numpy found successfully" FORCE)
    ENDIF ()

    IF (IS_DIRECTORY ${PYTHON_INCLUDE_DIR} AND EXISTS ${PYTHON_LIBRARIES})
      MESSAGE (STATUS "Python version: ${PYTHON_VERSION}")
      SET (PYTHON_FOUND TRUE)
    ELSE (IS_DIRECTORY ${PYTHON_INCLUDE_DIR} AND EXISTS ${PYTHON_LIBRARIES})
      MESSAGE (STATUS "Python not found in PYTHONHOME = $ENV{PYTHONHOME}")
      MESSAGE (STATUS "Searching in standard locations...")
      SET (PYTHON_FOUND FALSE)
    ENDIF (IS_DIRECTORY ${PYTHON_INCLUDE_DIR} AND EXISTS ${PYTHON_LIBRARIES})
    
  ENDIF (IS_DIRECTORY ${PYTHONROOT})
  
ENDIF (NOT PYTHON_FOUND)

IF (NOT PYTHON_FOUND)

  SET (PYTHON_FOUND TRUE)

  FIND_PACKAGE (PythonInterp QUIET)
  FIND_PACKAGE (PythonLibs QUIET) # secretly unsets PYTHON_FOUND

  IF (NOT PYTHON_EXECUTABLE)
    SET (PYTHON_FOUND FALSE)
  ELSE (NOT PYTHON_EXECUTABLE)
    SET (PYTHON_FOUND TRUE)
  ENDIF (NOT PYTHON_EXECUTABLE)
  
  EXECUTE_PROCESS (COMMAND ${PYTHON_EXECUTABLE} -V
    ERROR_VARIABLE PYTHON_VERSION_STRING
    ERROR_STRIP_TRAILING_WHITESPACE)
  SEPARATE_ARGUMENTS (PYTHON_VERSION_STRING)
  SET (PYTHON_VERSION "")
  FOREACH (item ${PYTHON_VERSION_STRING})
    SET (PYTHON_VERSION ${item})
  ENDFOREACH (item ${PYTHON_VERSION_STRING})

  # Search for numpy
  EXECUTE_PROCESS (COMMAND ${PYTHON_EXECUTABLE} -c "import numpy"
    RESULT_VARIABLE NUMPY_FOUND)
  IF (NUMPY_FOUND EQUAL 0)
    SET (NUMPY_FOUND TRUE)
  ELSE ()
    SET (NUMPY_FOUND FALSE)
  ENDIF ()

  IF (NUMPY_FOUND)
    SET (NUMPY_FOUND TRUE CACHE BOOL "Numpy found successfully" FORCE)
    EXECUTE_PROCESS (COMMAND ${PYTHON_EXECUTABLE} -c
      "import numpy; print numpy.get_include()"
      OUTPUT_VARIABLE NUMPY_INCLUDE_DIR
      OUTPUT_STRIP_TRAILING_WHITESPACE)
    SET (NUMPY_INCLUDE_DIR ${NUMPY_INCLUDE_DIR} CACHE STRING "Numpy directory")
  ELSE ()
    SET (NUMPY_FOUND FALSE CACHE BOOL "Numpy found successfully" FORCE)
  ENDIF ()

  MESSAGE (STATUS "Python version: ${PYTHON_VERSION}")

  SET (PYTHON_INCLUDE_DIR ${PYTHON_INCLUDE_PATH})
  IF (NOT EXISTS "${PYTHON_INCLUDE_DIR}/Python.h")
    MESSAGE (STATUS "Error: ${PYTHON_INCLUDE_DIR}/Python.h does not exist.\n")
    SET (PYTHON_FOUND FALSE)
    SET (PYTHON_CONFIG_ERROR TRUE)
  ENDIF (NOT EXISTS "${PYTHON_INCLUDE_DIR}/Python.h")

  SET (PYTHON_CPPFLAGS "-I${PYTHON_INCLUDE_DIR}")
  SET (PYTHON_LDFLAGS "${PYTHON_LIBRARIES}")

ENDIF (NOT PYTHON_FOUND)

IF (PYTHON_FOUND)
  MESSAGE (STATUS "  * binary:   ${PYTHON_EXECUTABLE}")
  MESSAGE (STATUS "  * includes: ${PYTHON_INCLUDE_DIR}")
  MESSAGE (STATUS "  * libs:     ${PYTHON_LIBRARIES}")
  IF (NUMPY_FOUND)
    MESSAGE (STATUS "  * numpy:    ${NUMPY_INCLUDE_DIR}")
  ENDIF ()
ENDIF (PYTHON_FOUND)
