# The name of our project is "FFTCW". CMakeLists files in this project can 
# refer to the root source directory of the project as ${FFTCW_SOURCE_DIR} and 
# to the root binary directory of the project as ${FFTCW_BINARY_DIR}. 
CMAKE_MINIMUM_REQUIRED(VERSION 2.6) 
PROJECT (FFTCW C) 

FIND_PATH(FFTW_INC_DIR
          NAMES QtGui
          PATHS /usr/include
                $ENV{HOME}/include
                /usr/local/include
                /opt/fftw/include
                /opt/fftw2/include
                /opt/fftw3/include
                /usr/local/fftw/include
                /usr/local/fftw2/include
                /usr/local/fftw3/include)


IF(NOT FFTW_INC_DIR)
   MESSAGE(FATAL_ERROR "Unable to find fftw include files")
ELSE(NOT FFTW_INC_DIR)
   MESSAGE(FATAL_ERROR "OK")
ENDIF(NOT FFTW_INC_DIR) 

ADD_DEFINITIONS( -Wall --pedantic -O3 )
ADD_EXECUTABLE(fftcw fftcw.c)


