# MRC_dump

A C++ and a python program that dump MRC image/movie/map slices into 8-bit grayscale PNG files. Very fast (compared to python/perl). Needs libpng.

The output PNG file is optimized for display. A linear scaling is applied with a default 1.5 sigma cutoff.

To compile:

g++ MapDump.cpp -std=c++11 -lpng -o MapDump.exe 

Example:

MapDump.exe *.mrcs -bin 3 -vbin 4 -sigma 2.5 -statframe 10
