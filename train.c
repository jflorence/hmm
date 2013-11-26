#include <stdlib.h>
#include "train.h"
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#define MAXCOUNTER 30


static inline void initalpha(double *alpha, int N, long T, int D)
{
	memset(alpha, 0, N*T);
	alpha[N*D] = 1.0;
	return;
}



static inline void initbeta(double *beta, int N, long T, int D)
{
	memset(beta, 0, N*T);
	for(int i = N*(T-D); i<(T-D)*(N+1);i++)
	{
		beta[i] = 1.0;
	}
}


static inline void compute_gamma_dist(double *dist, int length, double step, double shape, double scale)
{
	double x = 0;
	for(int i = 0; i<length; i++)
	{
		double G = tgamma(shape);
		dist[i] = (pow(x, shape-1)*exp(-x/scale))/(G*pow(scale, shape));
		x += step;
	}
}



void train(struct params *p, delays_mt *data, delay_mt ymax)
{

	const int N = p->N;
	const long T = data->length;
	//delay_mt *y = data->delay;


/*	double Ahat[N][N];
	double muhat[N];
	double sigmahat[N];
	double c[N];
*/	double alpha[N][T];
	double beta[N][T];
/*	double dshape[N];
	double dscale[N];
	double dshapehat[N];
	double dscalehat[N];
*/	int counter = 0;
	//double epsilon = 1.0;
	int D = p->D;
	assert(T>2*D);
	initalpha((double *) alpha, N, T, D);
	initbeta((double *) beta, N, T, D);

	double step = 0.001;
	double dstep = 1;
	long size_dist = ceil(ymax/step)+1;
	long size_dur_dist = D/dstep;
	double **distribs = malloc(sizeof(double *)*N);
	double **dur_dists = malloc(sizeof(double *)*N);
	if(distribs == NULL || dur_dists == NULL)
	{
		fprintf(stderr, "Couldn't malloc\n");
		exit(EXIT_FAILURE);
	}
	for(int i = 0; i< N; i++)
	{
		distribs[i] = malloc(sizeof(double)*size_dist);
		dur_dists[i] = malloc(sizeof(double)*size_dur_dist);
		if(distribs[i] == NULL || dur_dists[i] == NULL)
		{
			fprintf(stderr, "Couldn't malloc");
			exit(EXIT_FAILURE);
		}
	}


	printf("Computing gamma distributions...\n");
	for(int i=0; i<N; i++)
	{
		compute_gamma_dist(distribs[i], size_dist, step, p->shape[i], p->scale[i]);
		
		printf("shape: %f, scale: %f\n", p->shape[i], p->scale[i]);
		
		plot1D(distribs[i], size_dist);
	}


	//or likelihood...
//	while(counter < MAXCOUNTER)
	{
		/*compute_distributions();
		compute_duration_dist();
		compute_alphas();
		compute_betas();
		compute_A();
		compute_mu_sigma();
		//compute gamma params
		compute_duration_params();
		compute_duration_dist_params();

		//update A, mu, sigma...
		//update log likelihood
		*/
	}





}























































