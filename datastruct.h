#ifndef DATASTRUCT_H
#define DATASTRUCT_H

typedef double time_mt;
typedef double delay_mt;



struct delays_mt
{
	time_mt *time;
	delay_mt *delay;
	double length;
};
typedef struct delays_mt delays_mt;


struct params
{
	int N;
	double *A;
	double *shape;
	double *scale;
	double *mu;
	double *sigma;
	double *steady;
	double *pi;
	long D;

	//also, the params of the duration distributions:
	double *dmu;
	double *dsigma;
	double *dshape;
	double *dscale;
};












#endif
