#ifndef IO_H
#define IO_H

#include "datastruct.h"

#define VERBOSE

#define DEFAULTINPUT "delays.txt"
#define OUTPUTPARAMS "params.txt"



void parse(int argc, char *argv[]);
delays_mt getinput(void);
void dispstr(int level, char *str, ...);
void plot1D(float_mt *data, long long length);
void write_alpha_beta(float_mt *alpha, float_mt *beta, int N, long long L);
void write_results(struct params *p);
#endif

