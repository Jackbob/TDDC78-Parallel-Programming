#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <ctime>
#include "ppmio.h"
#include "thresfilter.h"
#include <mpi.h>
#include <iostream>

int main (int argc, char ** argv) {
    MPI_Init(nullptr, nullptr);
    int rank{}, world{}, root{0};

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world);

    int xsize, ysize, colmax;
    unsigned char* src, *newsrc;
    struct timespec stime{}, etime{};

    if(rank == root) {
        src = new unsigned char[MAX_PIXELS*3];
        newsrc = new unsigned char[MAX_PIXELS * 3];
        /* Take care of the arguments */

        if (argc != 3) {
            std::cerr << "Usage " << argv[0] << "infile outfile" << std::endl;
            return 0;
        }

        /* read file */
        if (read_ppm(argv[1], &xsize, &ysize, &colmax, (char *) src) != 0)
            return 0;

        if (colmax > 255) {
            std::cerr << "Too large maximum color-component value" << std::endl;
            return 0;
        }
    }

    std::cout << "Has read the image, calling filter" << std::endl;

    /*****************************************************
    * MPI PHASE                                          *
    *                                                    *
    ******************************************************/
    MPI_Bcast(&xsize, 1, MPI_INT, root, MPI_COMM_WORLD);
    MPI_Bcast(&ysize, 1, MPI_INT, root, MPI_COMM_WORLD);

    int ysplit = ysize/world;
    int rest = ysize%world;
    int give = rest==0 ? 0 : 1;
    int steal = world - 1 - rest;
    auto sendcounts = new int[world];
    auto displace = new int[world];
    for(int i=0; i<world; i++) {
        sendcounts[i] = (ysplit+give) * xsize * 3;
        displace[i] = (ysplit+give+1) * xsize * 3 * i;
    }

    if(rest != 0)
        sendcounts[world-1] = (ysplit - steal) * xsize * 3;

    auto dst = new unsigned char[sendcounts[rank]];

    MPI_Scatterv(src, sendcounts, displace, MPI_UNSIGNED_CHAR, dst, sendcounts[root], MPI_UNSIGNED_CHAR, root, MPI_COMM_WORLD);


    clock_gettime(CLOCK_REALTIME, &stime);

    //thresfilter(xsize, ysize, src);

    clock_gettime(CLOCK_REALTIME, &etime);
    MPI_Gather(dst, sendcounts[rank], MPI_UNSIGNED_CHAR, newsrc, sendcounts[root], MPI_UNSIGNED_CHAR, root, MPI_COMM_WORLD);

    std::cout << "Filtering took: " << (etime.tv_sec  - stime.tv_sec) +
                                       1e-9*(etime.tv_nsec  - stime.tv_nsec) << std::endl;

    /* write result */
    std::cout << "Writing output file" << std::endl;
    printf("Writing output file\n");

    if(write_ppm (argv[2], xsize, ysize, (char *)src) != 0)
      exit(1);

    if(rank == root){
        delete [] src;
    }
    delete [] sendcounts;
    delete [] displace;
    return(0);
}
