PROJECT(WCSLIB)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.0)

# Allow the developer to select if Dynamic or Static libraries are built
OPTION (BUILD_SHARED_LIBS "Build Shared Libraries" ON)
SET(CMAKE_MACOSX_RPATH OFF)
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}")
SET(LIB_DESTINATION "${CMAKE_INSTALL_PREFIX}/lib${LIB_SUFFIX}")
SET(INCLUDE_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/include/")
SET(INCLUDE_DIRS
    "${CMAKE_SOURCE_DIR}"
    "${CMAKE_SOURCE_DIR}/C"
)

# Define project version
SET(${PROJECT_NAME}_MAJOR_VERSION 5)
SET(${PROJECT_NAME}_MINOR_VERSION 15)
SET(${PROJECT_NAME}_VERSION ${${PROJECT_NAME}_MAJOR_VERSION}.${${PROJECT_NAME}_MINOR_VERSION})
SET(LIB_NAME wcslib)

# run ./configure for setting config file
if(EXISTS "${CMAKE_SOURCE_DIR}/wcsconfig.h")
  MESSAGE(STATUS "wcslib already configured...")
else()
  execute_process(COMMAND ./configure WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
endif()

#run lex files
if(EXISTS "${CMAKE_SOURCE_DIR}/C/fitshdr.c")
  MESSAGE(STATUS "flex fitshdr already created")
else()
  execute_process(COMMAND flex ${CMAKE_SOURCE_DIR}/C/fitshdr.l WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/C")
endif()

if(EXISTS "${CMAKE_SOURCE_DIR}/C/wcsbth.c")
  MESSAGE(STATUS "flex wcsbth already created")
else()
  execute_process(COMMAND flex ${CMAKE_SOURCE_DIR}/C/wcsbth.l WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/C")
endif()

if(EXISTS "${CMAKE_SOURCE_DIR}/C/wcspih.c")
  MESSAGE(STATUS "flex wcspih already created")
else()
  execute_process(COMMAND flex ${CMAKE_SOURCE_DIR}/C/wcspih.l WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/C")
endif()

if(EXISTS "${CMAKE_SOURCE_DIR}/C/wcsulex.c")
  MESSAGE(STATUS "flex wcsulex already created")
else()
  execute_process(COMMAND flex ${CMAKE_SOURCE_DIR}/C/wcsulex.l WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/C")
endif()

if(EXISTS "${CMAKE_SOURCE_DIR}/C/wcsutrn.c")
  MESSAGE(STATUS "flex wcsutrn already created")
else()
  execute_process(COMMAND flex ${CMAKE_SOURCE_DIR}/C/wcsutrn.l WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/C")
endif()

SET (LIB_TYPE STATIC)
IF (BUILD_SHARED_LIBS)
  SET (LIB_TYPE SHARED)
ENDIF (BUILD_SHARED_LIBS)

SET(H_FILES
    wcsconfig.h   wcsconfig_f77.h   wcsconfig_tests.h   wcsconfig_utils.h
    C/cel.h   C/dis.h
    C/fitshdr.h
    C/getwcstab.h
    C/lin.h   C/log.h
    C/prj.h   C/spc.h   C/sph.h   C/spx.h   C/tab.h   C/wcs.h
    C/wcserr.h   C/wcsfix.h   C/wcshdr.h   C/wcsprintf.h
    C/wcstrig.h
    C/wcsunits.h
    C/wcsutil.h
)

SET(SRC_FILES
    C/cel.c   C/dis.c
    C/fitshdr.c
    C/getwcstab.c
    C/lin.c   C/log.c
    C/prj.c   C/spc.c   C/sph.c   C/spx.c   C/tab.c   C/wcs.c
    C/wcsbth.c
    C/wcserr.c   C/wcsfix.c   C/wcshdr.c
    C/wcspih.c
    C/wcsprintf.c
    C/wcstrig.c
    C/wcsulex.c
    C/wcsulex.c
    C/wcsunits.c
    C/wcsutil.c
    C/wcsutrn.c
)

INCLUDE_DIRECTORIES(${INCLUDE_DIRS})

FIND_LIBRARY(CFITSIO_LIBRARY cfitsio
  PATHS ${CFITSIO_LIBRARY_DIR} NO_DEFAULT_PATH
  )

ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${H_FILES} ${SRC_FILES})
TARGET_LINK_LIBRARIES(${LIB_NAME} ${CFITSIO_LIBRARY})

SET_TARGET_PROPERTIES(${LIB_NAME} PROPERTIES VERSION ${${PROJECT_NAME}_VERSION} SOVERSION ${${PROJECT_NAME}_MAJOR_VERSION})
install(TARGETS ${LIB_NAME} LIBRARY DESTINATION ${LIB_DESTINATION})
install(FILES ${H_FILES} DESTINATION ${INCLUDE_INSTALL_DIR} COMPONENT Devel)

ENABLE_TESTING()
