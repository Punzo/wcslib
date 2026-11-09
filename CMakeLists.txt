# CMake build for WCSLIB (C library only).
#
# WCSLIB upstream builds with autotools; this file provides an equivalent
# CMake build of the C library so it can be used as a superbuild dependency
# (e.g. by SlicerAstro). Fortran wrappers, PGSBOX, utilities and tests are
# not built.
cmake_minimum_required(VERSION 3.20...3.31)

project(WCSLIB VERSION 8.10 LANGUAGES C)

include(CheckSymbolExists)
include(CheckTypeSize)
include(GNUInstallDirs)

option(BUILD_SHARED_LIBS "Build shared libraries" ON)
option(WCSLIB_USE_CFITSIO "Build getwcstab (requires CFITSIO)" ON)

set(WCSLIB_VERSION ${PROJECT_VERSION})
set(WCSLIB_LIBRARY_NAME "wcslib" CACHE STRING "Name of the built library")
mark_as_advanced(WCSLIB_LIBRARY_NAME)

if(NOT DEFINED CMAKE_C_STANDARD)
  # C11 is required for the thread-safe atomic spinlock in wcserr.c.
  set(CMAKE_C_STANDARD 11)
endif()
set(CMAKE_C_STANDARD_REQUIRED ON)

#-----------------------------------------------------------------------------
# Configuration checks (mirror configure.ac)

# 64-bit integer type
check_type_size("long long int" SIZEOF_LONG_LONG_INT)
check_type_size("long int" SIZEOF_LONG_INT)
if(SIZEOF_LONG_LONG_INT EQUAL 8)
  set(WCSLIB_INT64 "long long int")
elseif(SIZEOF_LONG_INT EQUAL 8)
  set(WCSLIB_INT64 "long int")
else()
  set(WCSLIB_INT64 "int")
endif()

# sincos()
set(CMAKE_REQUIRED_DEFINITIONS -D_GNU_SOURCE)
set(CMAKE_REQUIRED_LIBRARIES m)
check_symbol_exists(sincos "math.h" HAVE_SINCOS)
unset(CMAKE_REQUIRED_DEFINITIONS)
unset(CMAKE_REQUIRED_LIBRARIES)

# Thread-local storage specifier (used by wcsprintf.c)
if(MSVC)
  set(WCSLIB_TLS "__declspec(thread)")
else()
  set(WCSLIB_TLS "_Thread_local")
endif()

set(_staged_include_dir ${CMAKE_CURRENT_BINARY_DIR}/include)
configure_file(CMake/wcsconfig.h.in ${_staged_include_dir}/wcslib/wcsconfig.h @ONLY)

#-----------------------------------------------------------------------------
# Sources

set(_flexed fitshdr wcsbth wcspih wcsulex wcsutrn)
set(_modules cel dis lin log prj spc sph spx tab wcs wcserr wcsfix wcshdr
             wcsprintf wcstrig wcsunits wcsutil)

set(WCSLIB_SOURCES)
foreach(m IN LISTS _modules)
  list(APPEND WCSLIB_SOURCES C/${m}.c)
endforeach()
foreach(m IN LISTS _flexed)
  list(APPEND WCSLIB_SOURCES C/flexed/${m}.c)
endforeach()

set(WCSLIB_PUBLIC_HEADERS
  C/cel.h C/dis.h C/fitshdr.h C/lin.h C/log.h C/prj.h C/spc.h C/sph.h
  C/spx.h C/tab.h C/wcs.h C/wcserr.h C/wcsfix.h C/wcshdr.h C/wcslib.h
  C/wcsmath.h C/wcsprintf.h C/wcstrig.h C/wcsunits.h C/wcsutil.h
  C/wtbarr.h C/wcslimits.h)

if(WCSLIB_USE_CFITSIO)
  list(APPEND WCSLIB_SOURCES C/getwcstab.c)
  list(APPEND WCSLIB_PUBLIC_HEADERS C/getwcstab.h)

  find_path(CFITSIO_INCLUDE_DIR fitsio.h HINTS ${CFITSIO_INCLUDE_DIR})
  find_library(CFITSIO_LIBRARY cfitsio HINTS ${CFITSIO_LIBRARY_DIR})
  if(NOT CFITSIO_INCLUDE_DIR OR NOT CFITSIO_LIBRARY)
    message(FATAL_ERROR "WCSLIB_USE_CFITSIO is ON but CFITSIO was not found. "
                        "Set CFITSIO_INCLUDE_DIR and CFITSIO_LIBRARY_DIR.")
  endif()
endif()

# Stage public headers in the build tree as <wcslib/xxx.h>, so that
# dependent projects can consume the build tree directly.
foreach(h IN LISTS WCSLIB_PUBLIC_HEADERS)
  get_filename_component(_name ${h} NAME)
  configure_file(${h} ${_staged_include_dir}/wcslib/${_name} COPYONLY)
endforeach()

#-----------------------------------------------------------------------------
# Library

add_library(wcslib ${WCSLIB_SOURCES})
add_library(WCSLIB::wcslib ALIAS wcslib)
set_target_properties(wcslib PROPERTIES
  OUTPUT_NAME ${WCSLIB_LIBRARY_NAME}
  VERSION ${PROJECT_VERSION}
  SOVERSION ${PROJECT_VERSION_MAJOR}
  C_VISIBILITY_PRESET default
  POSITION_INDEPENDENT_CODE ON)

target_compile_definitions(wcslib PRIVATE "WCSLIB_TLS=${WCSLIB_TLS}")
if(HAVE_SINCOS)
  target_compile_definitions(wcslib PRIVATE _GNU_SOURCE)
endif()
if(WIN32)
  # Flex-generated scanners include <unistd.h> unless told otherwise.
  target_compile_definitions(wcslib PRIVATE YY_NO_UNISTD_H)
  # WCSLIB has no export annotations: export all functions from the DLL.
  # Note: global *data* (e.g. the *_errmsg arrays) is not exported this way;
  # consumers should use functions (e.g. wcserr/wcsperr) instead.
  set_target_properties(wcslib PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS ON)
endif()
if(MSVC)
  target_compile_definitions(wcslib PRIVATE _CRT_SECURE_NO_WARNINGS _USE_MATH_DEFINES)
  # wcserr.c uses <stdatomic.h> when available in C11 mode.
  if(MSVC_VERSION GREATER_EQUAL 1935)
    target_compile_options(wcslib PRIVATE /experimental:c11atomics)
  endif()
endif()

target_include_directories(wcslib
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/C
    ${_staged_include_dir}/wcslib
  PUBLIC
    $<BUILD_INTERFACE:${_staged_include_dir}>
    $<BUILD_INTERFACE:${_staged_include_dir}/wcslib>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/wcslib>)

if(WCSLIB_USE_CFITSIO)
  target_include_directories(wcslib PUBLIC $<BUILD_INTERFACE:${CFITSIO_INCLUDE_DIR}>)
  target_link_libraries(wcslib PUBLIC ${CFITSIO_LIBRARY})
endif()

if(NOT WIN32)
  target_link_libraries(wcslib PRIVATE m)
endif()

#-----------------------------------------------------------------------------
# Install

install(TARGETS wcslib EXPORT WCSLIBTargets
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT RuntimeLibraries
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT RuntimeLibraries
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development)

install(DIRECTORY ${_staged_include_dir}/wcslib
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} COMPONENT Development)

install(EXPORT WCSLIBTargets NAMESPACE WCSLIB::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/WCSLIB COMPONENT Development)

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/WCSLIBConfigVersion.cmake
  COMPATIBILITY SameMajorVersion)
file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/WCSLIBConfig.cmake
  "include(\"\${CMAKE_CURRENT_LIST_DIR}/WCSLIBTargets.cmake\")\n")
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/WCSLIBConfig.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/WCSLIBConfigVersion.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/WCSLIB COMPONENT Development)
