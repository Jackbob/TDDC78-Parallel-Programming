#include <mpi.h>
#include <iostream>
#include "blurfilter.h"
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstring>

#define MAX_X 1.33
#define Pi 3.14159

/* maximum number of pixels in a picture */
#define MAX_PIXELS (3000*3000)

int read_ppm (const char * fname, int * xpix, int * ypix, int * max, unsigned char * data);
int write_ppm (const char * fname, int xpix, int ypix, unsigned char * data);
void get_gauss_weights(int n, double* weights_out);

int main(int argc, char *argv[]) {

    MPI_Init(nullptr, nullptr);
    int rank{}, world{}, root{0};

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world);

    std::cout << "Rank " << rank << " out of " << world << "\n";
    int radius;
    int xsize, ysize, colmax;
    auto *src = new unsigned char[MAX_PIXELS * 3];

#define MAX_RAD 1000
    struct timespec stime{}, etime{};
    double w[MAX_RAD];


    if(rank == root) {


        /* Take care of the arguments */
        if (argc != 4) {
            fprintf(stderr, "Usage: %s radius infile outfile\n", argv[0]);
            exit(1);
        }
        radius = atoi(argv[1]);
        if ((radius > MAX_RAD) || (radius < 1)) {
            fprintf(stderr, "Radius (%d) must be greater than zero and less then %d\n", radius, MAX_RAD);
            exit(1);
        }

        /* read file */
        if (read_ppm(argv[2], &xsize, &ysize, &colmax, src) != 0)
            exit(1);

        if (colmax > 255) {
            fprintf(stderr, "Too large maximum color-component value\n");
            exit(1);
        }


        printf("Has read the image, generating coefficients\n");

        clock_gettime(CLOCK_REALTIME, &stime);
    }

    /* filter */
    get_gauss_weights(radius, w);

    MPI_Bcast(&xsize, 1, MPI_INT, root, MPI_COMM_WORLD);
    MPI_Bcast(&ysize, 1, MPI_INT, root, MPI_COMM_WORLD);
    MPI_Bcast(&radius, 1, MPI_INT, root, MPI_COMM_WORLD);
    MPI_Bcast(src, xsize*ysize*3, MPI_UNSIGNED_CHAR, root, MPI_COMM_WORLD);

    auto split = (int)std::ceil(ysize/world);
    auto* dst = new unsigned char[split*xsize*3];

    printf("Calling filter\n");

    blurfilter(xsize, ysize, src, dst, radius, w);


    MPI_Gather(dst, split*xsize*3, MPI_UNSIGNED_CHAR, src, split*xsize*3, MPI_UNSIGNED_CHAR, root, MPI_COMM_WORLD);

    if(rank == root) {
        /* write result */
        printf("Writing output file\n");

        clock_gettime(CLOCK_REALTIME, &etime);
        printf("Filtering took: %g secs\n", (etime.tv_sec - stime.tv_sec) +
                                            1e-9 * (etime.tv_nsec - stime.tv_nsec));

        if (write_ppm(argv[3], xsize, ysize, src) != 0)
            return 1;
    }

    delete(src);
    MPI_Finalize();
    return (0);
}






/* Generate an array of weights for the gaussian filter. */
/* Input: n - number of weights to generate              */
/* Output: weights_out - array of weights. The element [0] */
/*  should be used for the central pixel, elements [1..n] */
/*  should be used for the pixels on a distance [1..n]  */
void get_gauss_weights(int n, double* weights_out) {
    double x;
    int i;

    for( i = 0; i < n+1; i++) {
        x = (double)i * MAX_X/n;
        weights_out[i] = exp(-x*x * Pi);
    }
}


/* Function: read_ppm - reads data from an image file in PPM format.
   Input: fname - name of an image file in PPM format to read.
   Output:
      xpix, ypix - size of the image in x & y directions
      max - maximum intensity in the picture
      data - color data array. MUST BE PREALLOCATED to at least MAX_PIXELS*3 bytes.
   Returns: 0 on success.
 */
int read_ppm (const char * fname, int * xpix, int * ypix, int * max, unsigned char * data) {
    char ftype[40];
    char ctype[40] = "P6";
    char line[80];
    int err;

    FILE * fp;
    err = 0;

    if (fname == nullptr) fname = "\0";
    fp = fopen (fname, "r");
    if (fp == nullptr) {
        fprintf (stderr, "read_ppm failed to open %s: %s\n", fname,
                 strerror (err));
        return 1;
    }

    fgets(line, 80, fp);
    sscanf(line, "%s", ftype);

    while (fgets(line, 80, fp) && (line[0] == '#'));

    sscanf(line, "%d%d", xpix, ypix);
    fscanf(fp, "%d\n", max);

    if(*xpix * *ypix > MAX_PIXELS) {
        fprintf (stderr, "Image size is too big\n");
        return 4;
    };


    if (strncmp(ftype, ctype, 2) == 0) {
        if (fread (data, sizeof (char), *xpix * *ypix * 3, fp) !=
            *xpix * *ypix * 3) {
            perror ("Read failed");
            return 2;
        }
    } else {
        fprintf (stderr, "Wrong file format: %s\n", ftype);
    }

    if (fclose (fp) == EOF) {
        perror ("Close failed");
        return 3;
    }


    return 0;

}

/* Function: write_ppm - write out an image file in PPM format.
   Input:
      fname - name of an image file in PPM format to write.
      xpix, ypix - size of the image in x & y directions
      data - color data.
   Returns: 0 on success.
 */
int write_ppm (const char * fname, int xpix, int ypix, unsigned char * data) {

    FILE * fp;
    int err = 0;

    if (fname == nullptr) fname = "\0";
    fp = fopen (fname, "w");
    if (fp == nullptr) {
        fprintf (stderr, "write_ppm failed to open %s: %s\n", fname,
                 strerror (err));
        return 1;
    }

    fprintf (fp, "P6\n");
    fprintf (fp, "%d %d 255\n", xpix, ypix);
    if (fwrite (data, sizeof (char), xpix*ypix*3, fp) != xpix*ypix*3) {
        perror ("Write failed");
        return 2;
    }
    if (fclose (fp) == EOF) {
        perror ("Close failed");
        return 3;
    }
    return 0;
}

