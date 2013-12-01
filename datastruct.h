#ifndef DATASTRUCT_H
#define DATASTRUCT_H



#define DOUBLE  0
#define LDOUBLE 1
#define FLOAT   2

#define TYPET DOUBLE
#define TYPED LDOUBLE

#if TYPET == DOUBLE
typedef double time_mt;
#endif

#if TYPED == DOUBLE
typedef double delay_mt;
typedef double float_mt;
#elif TYPED == LDOUBLE
typedef long double delay_mt;
typedef long double float_mt;
#endif

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
