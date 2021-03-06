# The name of our project is "FFTCW". CMakeLists files in this project can 
# refer to the root source directory of the project as ${FFTCW_SOURCE_DIR} and 
# to the root binary directory of the project as ${FFTCW_BINARY_DIR}. 
CMAKE_MINIMUM_REQUIRED(VERSION 2.6) 
PROJECT (FFTCW C) 



SET(MYPATHS      /usr/include 
                $ENV{HOME}/include
                /usr/local/include
                /opt/fftw/include
                /opt/fftw3/include
                /usr/local/fftw/include
                /usr/local/fftw3/include)

FIND_PATH(FFTW_INC_DIR NAMES fftw3.h PATHS ${MYPATHS})


IF(NOT FFTW_INC_DIR)
   MESSAGE(FATAL_ERROR "Unable to find fftw include files")
ELSE(NOT FFTW_INC_DIR)
   MESSAGE(STATUS "Found fftw headers in ${FFTW_INC_DIR}")
   STRING(REGEX REPLACE "/fftw3.h" "" FFTW_INC_DIR "${FFTW_INC_DIR}")
   STRING(REGEX REPLACE "/include" "/lib" FFTW_LIB_DIR "${FFTW_INC_DIR}")
   FIND_PATH(FFTW_THREADS NAMES libfftw3_threads.a PATHS ${FFTW_LIB_DIR})
   IF(NOT FFTW_THREADS)
       SET(FFTW_LIBS fftw3)
       MESSAGE(STATUS "Unable to find threaded fftw version  in ${FFTW_THREADS}")
   ELSE(NOT FFTW_THREADS)
       SET(FFTW_LIBS fftw3 fftw3_threads)
       MESSAGE(STATUS "Found threaded fftw version in ${FFTW_THREADS}")
       ADD_DEFINITIONS( -D_USE_THREADS)

   ENDIF(NOT FFTW_THREADS) 
ENDIF(NOT FFTW_INC_DIR) 

INCLUDE_DIRECTORIES(${FFTW_INC_DIR})
LINK_DIRECTORIES(${FFTW_LIB_DIR})

ADD_DEFINITIONS( -Wall --pedantic -O3 )
ADD_EXECUTABLE(fftcw fftcw.c)
ADD_EXECUTABLE(fft fft.c)
TARGET_LINK_LIBRARIES(fftcw ${FFTW_LIBS})
TARGET_LINK_LIBRARIES(fft ${FFTW_LIBS})



