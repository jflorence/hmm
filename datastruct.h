#ifndef DATASTRUCT_H
#define DATASTRUCT_H

//typedef double  time_mt;
//typedef double delay_mt;
//typedef double float_mt;
#define time_mt  double
#define delay_mt double
#define float_mt double


struct delays_mt
{
	time_mt *time;
	delay_mt *delay;
	long long length;
};
typedef struct delays_mt delays_mt;


struct params
{
	int N;
	float_mt *A;
	float_mt *shape;
	float_mt *scale;
	float_mt *mu;
	float_mt *sigma;
	float_mt *steady;
	float_mt *pi;
	long D;

	//also, the params of the duration distributions:
	float_mt *dmu;
	float_mt *dsigma;
	float_mt *dshape;
	float_mt *dscale;
};












#endif
