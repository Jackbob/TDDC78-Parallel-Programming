cmake_minimum_required(VERSION 3.9)
project(Lab1)

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

SET(CMAKE_C_COMPILER mpicc)
SET(CMAKE_CXX_COMPILER mpicxx)

set(CMAKE_CXX_STANDARD 17)

add_executable(Lab1
        blur_fproxy.c
        blurfilter.c
        blurfilter.h
        blurmain.c
        gaussw.c
        gaussw.h
        ppm.c
        ppm.h
        ppmio.c
        ppmio.h
        thres_fproxy.c
        thresfilter.c
        thresfilter.h
        threshold.c
        thresmain.c
        )
