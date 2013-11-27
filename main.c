#include <stdlib.h>
#include <stdio.h>

#include "io.h"
#include "datastruct.h"
#include "aux.h"
#include "train.h"


#define MINDELAY 0.001
#define Ni 2
#define Di 2 //should be 4000
static void initparams(struct params *p, delay_mt ymax)
{
	p->N = Ni;
	p->A = malloc(sizeof(float_mt)*p->N);
	p->shape = malloc(sizeof(float_mt)*p->N);
	p->scale = malloc(sizeof(float_mt)*p->N);
	p->mu = malloc(sizeof(float_mt)*p->N);
	p->sigma = malloc(sizeof(float_mt)*p->N);
	p->dmu = malloc(sizeof(float_mt)*p->N);
	p->dsigma = malloc(sizeof(float_mt)*p->N);
	p->pi = malloc(sizeof(float_mt)*p->N);
	p->steady = malloc(sizeof(float_mt)*p->N);

	if(p->A == NULL || p->shape == NULL || p->scale == NULL || p->mu == NULL ||
		p->sigma == NULL || p->dmu == NULL || p->dsigma == NULL || 
		p->pi == NULL || p->steady == NULL)
	{
		fprintf(stderr, "Couldn't malloc\n");
		exit(EXIT_FAILURE);
	}
	
	for(int i=0;i<p->N; i++)
	{
		p->A[i] = 1/(p->N-1);
	}
	//initialize the diagonal to 0.
	for(int i=0; i<p->N; i++)
	{
		p->A[i] = 0;
		i = i+p->N+1;
	}

	for(int i=0;i<p->N; i++)
	{
		p->shape[i] = (i+1)*ymax/(p->N+1);
		p->scale[i] = 1.0;
		p->mu[i] = 0;
		p->dmu[i] = i+1;
		p->dsigma[i] = 1;
	}
	
	p->D = Di;

	
}


int main(int argc, char *argv[])
{
	printf("Welcome to the HMM program!\ni");
	parse(argc, argv);
	printf("Loading file...\n");
	delays_mt input = getinput();
	printf("Sorting inputs...\n");
	mysort(&input);

	//this offsets the timestamps to zero, and finds the min and max value of the delays.
	delay_mt ymin = input.delay[0];
	delay_mt ymax = ymin;
	delay_mt current;
	time_mt tmin = input.time[0];
	for(int i = 0; i<input.length; i++)
	{
		input.time[i] = input.time[i] - tmin;
		current = input.delay[i];
		if(current < ymin)
			ymin = current;
		if(current > ymax)
			ymax = current;
	}
	//now, we offset the delays to zero
	for(int i = 0;i<input.length;i++)
	{   //MINDELAY is an epsilon in order not to bug the algo
		input.delay[i] = input.delay[i] - ymin + MINDELAY;
	}


	struct params p;
	printf("Initializing parameters...\n");
	initparams(&p, ymax);
	printf("Training model...\n");
	train(&p, &input, ymax);



	return EXIT_SUCCESS;
}
