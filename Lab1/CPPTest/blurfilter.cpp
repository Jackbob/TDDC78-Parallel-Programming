/*
  File: blurfilter.c

  Implementation of blurfilter function.

 */
#include <cstdio>
#include "blurfilter.h"
#include "ppmio.h"
#include <mpi.h>
#include <ctgmath>


pixel* pix(pixel* image, const int xx, const int yy, const int xsize)
{
  int off = xsize*yy + xx;

#ifdef DBG
  if(off >= MAX_PIXELS) {
    fprintf(stderr, "\n Terribly wrong: %d %d %d\n",xx,yy,xsize);
  }
#endif
  return (image + off);
}

void blurfilter(const int xsize, const int ysize, unsigned char *src, unsigned char *dst, const int radius, const double *w) {

    int rank{}, world{}, root{0};
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world);

    auto split = (int)std::ceil(ysize/world);

    int from = split * rank;
    int to = from + split;

    int x{}, y{};

    for(int r = from; r<to; r++) {
        for(int c = 0; c<xsize; c++) {
            dst[xsize*(r-from)*3 + c*3 + 0] = 0;
            dst[xsize*(r-from)*3 + c*3 + 1] = 0;
            dst[xsize*(r-from)*3 + c*3 + 2] = 0;
            int n = 0;
            for (x = -radius; x <= radius; x++) {
                for (y = -radius; y <= radius; y++) {
                    if(!(c+x<0 || c+x>xsize || r+y<0 || r+y>ysize)){
                        dst[xsize*r*3 + c*3 + 0] += src[xsize*3*(from+r+y) + (c+x)*3 + 0];
                        dst[xsize*r*3 + c*3 + 1] += src[xsize*3*(from+r+y) + (c+x)*3 + 1];
                        dst[xsize*r*3 + c*3 + 2] += src[xsize*3*(from+r+y) + (c+x)*3 + 2];
                        n++;
                    }
                }
            }

            dst[xsize*r*3 + c*3 + 0] = dst[xsize*3*(from+r+y) + (c+x)*3 + 0] / n;
            dst[xsize*r*3 + c*3 + 1] = dst[xsize*3*(from+r+y) + (c+x)*3 + 1] / n;
            dst[xsize*r*3 + c*3 + 2] = dst[xsize*3*(from+r+y) + (c+x)*3 + 2] / n;

        }
    }

}
