# Original source is Copyright Olivier Parcollet 2014.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

# Modified and Extended by Florian Goth 2018.
#
# This module looks for the nfft:
# https://www-user.tu-chemnitz.de/~potts/nfft/
# It sets up the following variables:
#
# NFFT_FOUND               - True if some version of NFFT was found.
# NFFT_INCLUDE_DIR         - Path to the include files
# NFFT_LIBRARIES           - Path to the libraries
# NFFT_VERSION             - the version of NFFT found as a string
# NFFT_VERSION_MAJOR       - the major version number of NFFT
# NFFT_VERSION_MINOR       - the minor version number of NFFT
# NFFT_VERSION_DECIMAL     - e.g. version 1.6.1 is 10601
# Version checks can be done directly in the FindNFFT call.

SET(TRIAL_PATHS
 $ENV{NFFT_ROOT}/include
 ${NFFT_ROOT}/include
 ~/opt/include
 ~/.local/include
 ~/include
 /usr/include
 /usr/local/include
 /opt/local/include
 /opt/nfft/include
 /sw/include
 )
FIND_PATH(NFFT_INCLUDE_DIR nfft3.h ${TRIAL_PATHS} DOC "Include for NFFT")

SET(TRIAL_LIBRARY_PATHS
 $ENV{NFFT_ROOT}/lib
 ${NFFT_ROOT}/lib
 ~/opt/lib
 ~/.local/lib
 ~/lib
 /usr/lib 
 /usr/local/lib
 /opt/local/lib
 /opt/nfft/lib
 /sw/lib
 )

SET(NFFT_LIBRARIES "NFFT_LIBRARIES-NOTFOUND" CACHE STRING "NFFT library")

# Try to detect the lib
FIND_LIBRARY(NFFT_LIBRARIES nfft3 ${TRIAL_LIBRARY_PATHS} DOC "NFFT library")

mark_as_advanced(NFFT_INCLUDE_DIR)
mark_as_advanced(NFFT_LIBRARIES)

if(NOT (NFFT_LIBRARIES STREQUAL "NFFT_LIBRARIES-NOTFOUND"))
    SET(WORK_DIR ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/FindNFFT)
    SET(NFFT_VERSION_SRC 
    "#include <stdio.h>
    #include <nfft3.h>
    int main(){
    unsigned major, minor, patch;
    nfft_get_version(&major, &minor, &patch);
    printf(\"%u.%u.%u\", major, minor, patch);
    return 0;}")
    file(WRITE "${WORK_DIR}/nfft3test.c" "${NFFT_VERSION_SRC}")
    try_run(RUN_RESULT_VAR COMPILE_RESULT_VAR ${CMAKE_BINARY_DIR} "${WORK_DIR}/nfft3test.c"
    LINK_LIBRARIES ${NFFT_LIBRARIES}
    COMPILE_OUTPUT_VARIABLE comp
    RUN_OUTPUT_VARIABLE NFFT_VERSION)
    if(COMPILE_RESULT_VAR)
        if(NOT RUN_RESULT_VAR) # cmake interprets 1 as true. 
        # Get the major and minor version numbers
        string(REGEX REPLACE "\\." ";" NFFT_VERSION_LIST ${NFFT_VERSION})
        list(GET NFFT_VERSION_LIST 0 NFFT_VERSION_MAJOR)
        list(GET NFFT_VERSION_LIST 1 NFFT_VERSION_MINOR)
        list(GET NFFT_VERSION_LIST 2 NFFT_VERSION_PATCH)
#        string(REGEX MATCH "[0-9]*" NUMPY_VERSION_PATCH ${NUMPY_VERSION_PATCH})        
        else()
        message(FATAL_ERROR "[NFFT] Compilation successful, but could not run executable!")
        endif()
    else()
    # Compilation failed This happens for NFFT before 3.3.0 since it lacks the get_version call
    # Since the library was found we return 1.0.0 as version string.
    # For the future: in case precise version numbers would be required for NFFT < 3.3.0, we could ask packageconfig
    message(STATUS "[NFFT] Library found, but could not compile test! Assuming old version of NFFT")
    SET(NFFT_VERSION_MAJOR "1")
    SET(NFFT_VERSION_MINOR "0")
    SET(NFFT_VERSION_PATCH "0")
    SET(NFFT_VERSION "${NFFT_VERSION_MAJOR}.${NFFT_VERSION_MINOR}.${NFFT_VERSION_PATCH}")
    endif()
    math(EXPR NFFT_VERSION_DECIMAL
        "(${NFFT_VERSION_MAJOR} * 10000) + (${NFFT_VERSION_MINOR} * 100) + ${NFFT_VERSION_PATCH}")
endif()

FIND_PACKAGE_HANDLE_STANDARD_ARGS(NFFT
                                  REQUIRED_VARS NFFT_LIBRARIES NFFT_INCLUDE_DIR NFFT_VERSION
                                  VERSION_VAR NFFT_VERSION)
