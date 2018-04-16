cmake_minimum_required(VERSION 3.9)
project(Lab1)

set(CMAKE_CXX_STANDARD 11)

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
        thresmain.c)
