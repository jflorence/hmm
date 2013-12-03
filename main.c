#include <stdlib.h>
#include <stdio.h>
#include "io.h"
#include "datastruct.h"
#include "aux.h"
#include "train.h"

#define MINDELAY 0.001
#define Ni 3
#define Di 2//should be 4000

static void initparams(struct params *p, delay_mt ymax);
static void freeparams(struct params *p);
int main(int argc, char *argv[])
{
	dispstr(0,"Welcome to the HMM program!\n");

#ifndef __STDC_IEC_559__
	printf("WARNING: ISO/IEC 60559 not respected!\n");
#endif

	parse(argc, argv);
	dispstr(0,"Loading file...\n");
	delays_mt input = getinput();
	dispstr(0,"Sorting inputs...\n");
	mysort(&input);


	//this offsets the timestamps to zero, and finds the min and max value of the delays.
	delay_mt ymin = input.delay[0];
	delay_mt ymax = ymin;
	delay_mt current;
	time_mt tmin = input.time[0];
	for(long long i = 0; i<input.length; i++)
	{
		input.time[i] = input.time[i] - tmin;
		current = input.delay[i];
		if(current < ymin)
			ymin = current;
		if(current > ymax)
			ymax = current;
	}
	//now, we offset the delays to zero
	for(long long i = 0;i<input.length;i++)
	{   //MINDELAY is an epsilon in order not to bug the algo
		input.delay[i] = input.delay[i] - ymin + MINDELAY;
	}


	struct params p;
	dispstr(0,"Initializing Markov parameters...\n");
	initparams(&p, ymax);
	

	input.length = input.length/300;

	dispstr(0,"Training model...\n");
	train(&p, &input, ymax);

	dispstr(0,"Writing results...\n");
	write_results(&p);

	freeparams(&p);

	dispstr(0,"Done.\n");


	return EXIT_SUCCESS;
}

static void freeparams(struct params *p)
{
	free(p->A);
	free(p->shape);
	free(p->scale);
	free(p->mu);
	free(p->sigma);
	free(p->dmu);
	free(p->dsigma);
	free(p->dshape);
	free(p->scale);
	free(p->pi);
	free(p->steady);
}

static void initparams(struct params *p, delay_mt ymax)
{
	p->N = Ni;
	p->D = Di;
	p->A = malloc(sizeof(float_mt)*p->N*p->N);
	p->shape = malloc(sizeof(float_mt)*p->N);
	p->scale = malloc(sizeof(float_mt)*p->N);
	p->mu = malloc(sizeof(float_mt)*p->N);
	p->sigma = malloc(sizeof(float_mt)*p->N);
	p->dmu = malloc(sizeof(float_mt)*p->N);
	p->dsigma = malloc(sizeof(float_mt)*p->N);
	p->dshape = malloc(sizeof(float_mt)*p->N);
	p->dscale = malloc(sizeof(float_mt)*p->N);
	p->pi = malloc(sizeof(float_mt)*p->N);
	p->steady = malloc(sizeof(float_mt)*p->N);

	if(p->A == NULL || p->shape == NULL || p->scale == NULL || p->mu == NULL ||
		p->sigma == NULL || p->dmu == NULL || p->dsigma == NULL || 
		p->pi == NULL || p->steady == NULL)
	{
		fprintf(stderr, "Couldn't malloc\n");
		exit(EXIT_FAILURE);
	}

	if(p->N>1)
	{
		for(int i=0;i<p->N*p->N; i++)
		{
			p->A[i] = 1.0/(p->N-1);
		}
		//initialize the diagonal to 0.
		for(int i=0; i<p->N*p->N; i++)
		{
			p->A[i] = 0.0;
			i = i+p->N;
		}
	}
	else
	{
		p->A[0] = 1.0;
	}

	for(int i=0;i<p->N; i++)
	{
		//p->shape[i] = (i+1.0)*ymax/(p->N+1);
		p->shape[i] = i+1.0;
		p->scale[i] = 1.0;
		p->mu[i] = 0.0;
		p->dmu[i] = i+1.0;
		p->dsigma[i] = 1.0;
		p->dshape[i] = i+1.0;   //how to initialize those ?
		p->dscale[i] = i+1.0;
	}
	
}











