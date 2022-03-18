# -*- mode: cmake; cmake-tab-width: 4-*- # vim: et:

# build RPM or DEB packages of photopsline. for macOS use Homebrew

if(CMAKE_SYSTEM_NAME MATCHES "Linux")

  file(STRINGS /etc/os-release OSR
    REGEX ^ID_LIKE=)

  set(CPACK_PACKAGE_VERSION ${CMAKE_PROJECT_VERSION})
  set(CPACK_PACKAGE_NAME "${CMAKE_PROJECT_NAME}")
  set(CPACK_PACKAGE_RELEASE 1)
  set(CPACK_PACKAGE_CONTACT "https://github.com/icecube/photospline")
  set(CPACK_PACKAGE_VENDOR "IceCube SPNO")
  set(CPACK_PACKAGE_CHECKSUM SHA1)
  set(CPACK_PACKAGE_DESCRIPTION "photopsline")
  set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}-${CPACK_PACKAGE_RELEASE}.${CMAKE_SYSTEM_PROCESSOR}")

  set(CPACK_PACKAGING_INSTALL_PREFIX "/usr")

  if(OSR MATCHES "rhel")
    set(CPACK_GENERATOR RPM)
    set(CPACK_RPM_PACKAGE_LICENSE "BSD")
    set(CPACK_RPM_PACKAGE_REQUIRES "python >= 3, cfitsio, gsl, openblas, suitesparse")
    set(CPACK_RPM_BUILDREQUIRES "python-devel >= 3, cfitsio-devel, gsl-devel, openblas-devel, suitesparse-devel")
  elseif(OSR MATCHES "debian")
    set(CPACK_GENERATOR DEB)
    set(CPACK_DEBIAN_PACKAGE_DEPENDS "python3, libcfitsio9, libgsl27, libopenblas0, libamd2, libccolamd2, libcholmod3, libcolamd2, libspqr2")
  else()
    # unknown/unsupported package type
    return()
  endif()

  include(CPack)

endif()
