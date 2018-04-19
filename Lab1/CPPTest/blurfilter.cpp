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

void blurfilter(const int xsize, const int ysize, unsigned char* overlap_top, unsigned char* overlap_bot, unsigned char *dst, const int radius, const double *w) {

    int rank{}, world{}, root{0};
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world);

    auto split = (int)std::ceil(ysize/world);

    int x{}, y{};
    double wc{}, n{};


    for(int row = 0; row<split; row++) {
        for(int c = 0; c<xsize; c++) {
            unsigned char r = 0;
            unsigned char g = 0;
            unsigned char b = 0;

            n = w[0];
            for (x = -radius; x <= radius; x++) {
                for (y = -radius; y <= radius; y++) {
                    wc = w[std::max(abs(x), abs(y))];
                    if(!(c+x<0 || c+x>xsize)){
                        if( row+y < 0  && rank != root){
                            r += wc * overlap_top[xsize * 3*(radius+y) + (c+x)*3 + 0];
                            g += wc * overlap_top[xsize * 3*(radius+y) + (c+x)*3 + 1];
                            b += wc * overlap_top[xsize * 3*(radius+y) + (c+x)*3 + 2];
                            n += wc;
                        }
                        else if( row+y > ysize && rank != world-1){
                            r += wc * overlap_bot[xsize * 3*(radius+y) + (c+x)*3 + 0];
                            g += wc * overlap_bot[xsize * 3*(radius+y) + (c+x)*3 + 1];
                            b += wc * overlap_bot[xsize * 3*(radius+y) + (c+x)*3 + 2];
                            n += wc;
                        }
                        else {
                            r += wc * dst[xsize * 3*(row+y) + (c+x)*3 + 0];
                            g += wc * dst[xsize * 3*(row+y) + (c+x)*3 + 1];
                            b += wc * dst[xsize * 3*(row+y) + (c+x)*3 + 2];
                            n += wc;
                        }
                    }
                }
            }

            dst[xsize*row*3 + c*3 + 0] = static_cast<unsigned char>(r / n);
            dst[xsize*row*3 + c*3 + 1] = static_cast<unsigned char>(g / n);
            dst[xsize*row*3 + c*3 + 2] = static_cast<unsigned char>(b / n);

        }
    }


}
