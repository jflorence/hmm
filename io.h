#ifndef IO_H
#define IO_H

#define DEFAULTINPUT "delays.txt"
#include "datastruct.h"

void parse(int argc, char *argv[]);
delays_mt getinput(void);

void plot1D(double *data, long length);

#endif

