#include <iostream>
#include "thresfilter.h"
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <algorithm>
#include <array>
#include <iterator>
#include <pthread.h>
#include <thread>
#include <vector>

#define MAX_X 1.33
#define Pi 3.14159

/* maximum number of pixels in a picture */
#define MAX_PIXELS (3000*3000)

int read_ppm (const char * fname, int * xpix, int * ypix, int * max, unsigned char * data);
int write_ppm (const char * fname, int xpix, int ypix, unsigned char * data);

int main(int argc, char *argv[]) {
    unsigned int num_processor = std::thread::hardware_concurrency();
    std::cout << num_processor << std::endl;
    std::vector<pthread_t> threads{num_processor};
    
    int xsize, ysize, colmax;
    unsigned char* src;
    unsigned char* newsrc;

    struct timespec stime{}, etime{};

    src = new unsigned char[MAX_PIXELS * 3];
    newsrc = new unsigned char[MAX_PIXELS * 3];

    if (argc != 3) {
      std::cerr << "Usage: %s infile outfile\n" << argv[0] << std::endl;
      exit(1);
      }


      if(read_ppm (argv[1], &xsize, &ysize, &colmax, src) != 0)
          exit(1);

      if (colmax > 255) {
          std::cerr << "Too large maximum color-component value\n" << std::endl;
    }

    std::cout << "Has read the image, generating coefficients <<" std::endl;

    clock_gettime(CLOCK_REALTIME, &stime);

    int ysplit = ysize/num_processor;
    int rest = ysize%num_processor;
    std::vector<unsigned int> splitcounts(num_processor, ysplit);
    splitcounts[num_processor-1] += rest;
    std::copy(splitcounts.begin(), splitcounts.end(), std::ostream_iterator<unsigned int>(std::cout, " "));

    std::cout << "Calling filter" << std::endl;

    /*unsigned int local_sum{0},sum{0};

    for(int i = 0; i < splitcounts; i++){
      local_sum += src[i];
    }*/

    //unsigned char mean = static_cast<unsigned char>(sum/(xsize*ysize*3));

    //thresfilter(xsize, sendcounts[rank]/(xsize*3), src, newsrc, mean);

    std::cout << "Calling filter" << std::endl;

    clock_gettime(CLOCK_REALTIME, &etime);

    std::cout << "Filtering took" << (etime.tv_sec - stime.tv_sec) +
                                     1e-9 * (etime.tv_nsec - stime.tv_nsec) << " secs" << std::endl;

    if (write_ppm(argv[2], xsize, ysize, newsrc) != 0)
        return 1;

    delete[] src;
    delete[] newsrc;

    return (0);
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
