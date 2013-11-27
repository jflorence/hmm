#include <stdlib.h>
#include "train.h"
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "aux.h"


#define MAXCOUNTER 30



static int N;
static long T;
static long D;
static double step = 0.001;
static double dstep = 1;
static double *alpha;
static double *beta;
static double **ydensity;
static double **distribs;
static double **dur_dists;
static delay_mt *y;
static double *A;
static double *Ahat;
double *muhat;
double *sigmahat;
double *c;
double *dshape;
double *dscale;
double *dshapehat;
double *dscalehat;





static inline void initalpha(void);
static inline void initbeta(void);
static inline void compute_gamma_dist(double *dist, int length, double step, double shape, double scale);

static inline void density_of_y();
static inline void compute_alphas(void);
static inline void compute_betas(void);
static inline void compute_A(void);
static inline void compute_mu_sigma(void);
//compute gamma params
static inline void compute_duration_params(void);
static inline void compute_duration_dist_params(void);

void train(struct params *p, delays_mt *data, delay_mt ymax)
{

	N = p->N;
	T = data->length;
	y = data->delay;

	alpha = malloc(N*T*sizeof(double));
	beta = malloc(N*T*sizeof(double));
	if(alpha == NULL || beta == NULL)
	{
		fprintf(stderr, "Couldn't malloc\n");
		exit(EXIT_FAILURE);
	}

	int counter = 0;
	//double epsilon = 1.0;
	D = p->D;
	assert(T>2*D);
	initalpha();
	initbeta();

	long size_dist = ceil(ymax/step)+1;
	long size_dur_dist = D/dstep;
	distribs = malloc(sizeof(double *)*N);
	dur_dists = malloc(sizeof(double *)*N);
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

	ydensity = malloc(N*sizeof(double));
	if(ydensity == NULL)
	{
		fprintf(stderr, "Couldnt malloc\n");
		exit(EXIT_FAILURE);
	}
	for(int i =0; i<N; i++)
	{
		ydensity[i] = malloc(sizeof(double)*T);
		if(ydensity[i] == NULL)
		{
			fprintf(stderr, "Couldn't malloc\n");
			exit(EXIT_FAILURE);
		}
	}

	A = malloc(N*N*sizeof(double));
	memcpy(A, p->A, N*N*sizeof(double));
	c = malloc(T*sizeof(double));
	for(int i = 0; i<T; i++)
		c[i] = 1.0;
	Ahat = malloc(N*N*sizeof(double));
	muhat = malloc(sizeof(double)*N);
	sigmahat = malloc(sizeof(double)*N);
	dshape = malloc(sizeof(double)*N);
	dscale = malloc(sizeof(double)*N); 
	dshapehat = malloc(sizeof(double)*N); 
	dscalehat = malloc(sizeof(double)*N); 
	if(A == NULL || Ahat == NULL || muhat == NULL || sigmahat == NULL ||
		dshape == NULL || dscale == NULL || dshapehat == NULL || dscalehat == NULL)
	{
		fprintf(stderr, "Couldn't malloc\n");
		exit(EXIT_FAILURE);
	}
	
	
	//or likelihood...
//	while(counter < MAXCOUNTER)
	{
		printf("Starting iteration %d\n", ++counter);
		printf("Computing gamma distributions...\n");
		for(int i=0; i<N; i++)
		{
			compute_gamma_dist(distribs[i], size_dist, step, p->shape[i], p->scale[i]);
			compute_gamma_dist(dur_dists[i], size_dur_dist, dstep, p->dshape[i], p->dscale[i]);	
		
		//	plot1D(distribs[i], size_dist);
		}

		density_of_y();
		initalpha();
		compute_alphas();
		initbeta();
		compute_betas();
		compute_A();
		compute_mu_sigma();
		//compute gamma params
		compute_duration_params();
		compute_duration_dist_params();

		//update A, mu, sigma...
		//update log likelihood
	}


}

//alpha and beta should be stored so that the values for the same t are contiguous.
static inline void initalpha(void)
{
	memset(alpha, 0, N*T);
	alpha[N*D] = 1.0;
	return;
}



static inline void initbeta(void)
{
	memset(beta, 0, N*T);
	for(int i = N*(T-D); i<(T-D+1)*(N);i++)
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
	return;
}



void density_of_y(void)
{
	//we need the proba density for each sample of the y vector, according to the current distributions
	int index;
	for(int j=0; j<N; j++)
	{
		for(int k = 0; k<T; k++)
		{
			index = round(y[k]/step)+1; //this is an index on the vector of the pdf
			(ydensity[j])[k] = (distribs[j])[index];
		}
	}
}


static inline void compute_alphas(void)
{
	//Is this the best order for the loops ?
	register double term;
	
	for(int t = D; t<T-D; t++)
	{
		for(int j = 0; j<N; j++)
		{
			for(int d = 0; d<D; d++)
			{
				for(int i=0; i<N; i++)
				{
					term = alpha[(t-d-1)*N+i]*A(i*N+j)*(dur_dists[j])[d]*prod(&c[t-d],d-1);
					term *= prod(&(ydensity[j])[t-d],d);
					alpha[t*N+j] += term;
				}
			}
		}
		c[t] = 1/sum(&alpha[T*N], N);
		for(i=0; i<N; i++)
		{
			alpha[t*N + i] *= c[t];
		}
	}
}
static inline void compute_betas(void)
{
	register double term;
	long index;
	for(int t = T-D-1; t>=D; t--)
	{
		for(int i =0; i<N; i++)
		{
			index = t*N+i;
			for(int d=0; d<D; d++)
			{
				for(int j = 0; j<N; j++)
				{
					term = A[i*N+j]*(dur_dist[j])[d]*prod(c[t+1], d-1)*beta((t+d-1)*N+j);
					term *= prod(&(ydensity[j])[t+1],d);
					beta[index] += term;
				}
			}
			beta[index] *= c[t];
		}
	}
}
static inline void compute_A(void)
{
	double num;
	double term;
	double den;
	double *temp = malloc(sizeof(double)*T);
	for(int i=0;i<N;i++)
	{
		for(int j=0; j<N; j++)
		{
			num = 0;
			for(int t = D; t<T-D; t++)
			{
				for(int d=0; d<D; d++)
				{
					term = A[i*N+j]*(dur_dist[j])[d]*beta[(t+d)*N+j]*prod(c[t+1],d-1);
					term *= alpha[t*N+i]*prod(&(ydensity[j])[t+1],d);
					num += term;
				}
			temp[t]=alpha[t*N+i]*beta[t*N+i]/c[t];
			}
			dem = sum(&temp[D], T-2*D);
			Ahat[i*N+j] = num/den;
		}
	}
	free(temp);
}
static inline void compute_mu_sigma(void)
{
}
static inline void compute_duration_params(void)
{
}
static inline void compute_duration_dist_params(void)
{
}










