# MRC_dump

A C++ program that dumps MRC image/movie/map slices into 8-bit grayscale PNG files. Very fast (compared to python/perl). Needs libpng.

+/- 6-RMS cutoff is applied to filter out X-ray/flare pixels. The RMS, mean, max and min values are always recalculated from slice 0 and is used for calculating the scaling factor and cutoff values, which are applied to all slices.

To compile:

g++ MapDump.cpp -std=c++11 -lpng -o MapDump.exe 
