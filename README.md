# MRC_dump

Dumps MRC image/movie/map slices into 8-bit grayscale PNG files. Very fast. Needs libpng.

The result PNG show rescaled(linearly) images. 6-RMS cutoff is applied to filter out X-ray/flare pixels.

To compile:

g++ MapDump.cpp -std=c++11 -lpng -o MapDump.exe 
