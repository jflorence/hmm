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
static float_mt step = 0.001;
static float_mt dstep = 1;
static float_mt *alpha;
static float_mt *beta;
static float_mt **ydensity;
static float_mt **distribs;
static float_mt **dur_dists;
static delay_mt *y;
static float_mt *A;
static float_mt *Ahat;
float_mt *muhat;
float_mt *sigmahat;
float_mt *c;
float_mt *dshape;
float_mt *dscale;
float_mt *dshapehat;
float_mt *dscalehat;





static inline void initalpha(void);
static inline void initbeta(void);
static inline void compute_gamma_dist(float_mt *dist, int length, float_mt step, float_mt shape, float_mt scale);

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

	alpha = malloc(N*T*sizeof(float_mt));
	beta = malloc(N*T*sizeof(float_mt));
	if(alpha == NULL || beta == NULL)
	{
		fprintf(stderr, "Couldn't malloc\n");
		exit(EXIT_FAILURE);
	}

	int counter = 0;
	//float_mt epsilon = 1.0;
	D = p->D;
	assert(T>2*D);
	initalpha();
	initbeta();

	long size_dist = ceil(ymax/step)+1;
	long size_dur_dist = D/dstep;
	distribs = malloc(sizeof(float_mt *)*N);
	dur_dists = malloc(sizeof(float_mt *)*N);
	if(distribs == NULL || dur_dists == NULL)
	{
		fprintf(stderr, "Couldn't malloc\n");
		exit(EXIT_FAILURE);
	}
	for(int i = 0; i< N; i++)
	{
		distribs[i] = malloc(sizeof(float_mt)*size_dist);
		dur_dists[i] = malloc(sizeof(float_mt)*size_dur_dist);
		if(distribs[i] == NULL || dur_dists[i] == NULL)
		{
			fprintf(stderr, "Couldn't malloc");
			exit(EXIT_FAILURE);
		}
	}

	ydensity = malloc(N*sizeof(float_mt));
	if(ydensity == NULL)
	{
		fprintf(stderr, "Couldnt malloc\n");
		exit(EXIT_FAILURE);
	}
	for(int i =0; i<N; i++)
	{
		ydensity[i] = malloc(sizeof(float_mt)*T);
		if(ydensity[i] == NULL)
		{
			fprintf(stderr, "Couldn't malloc\n");
			exit(EXIT_FAILURE);
		}
	}

	A = malloc(N*N*sizeof(float_mt));
	memcpy(A, p->A, N*N*sizeof(float_mt));
	c = malloc(T*sizeof(float_mt));
	for(int i = 0; i<T; i++)
		c[i] = 1.0;
	Ahat = malloc(N*N*sizeof(float_mt));
	muhat = malloc(sizeof(float_mt)*N);
	sigmahat = malloc(sizeof(float_mt)*N);
	dshape = malloc(sizeof(float_mt)*N);
	dscale = malloc(sizeof(float_mt)*N); 
	dshapehat = malloc(sizeof(float_mt)*N); 
	dscalehat = malloc(sizeof(float_mt)*N); 
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


static inline void compute_gamma_dist(float_mt *dist, int length, float_mt step, float_mt shape, float_mt scale)
{
	float_mt x = 0;
	for(int i = 0; i<length; i++)
	{
		float_mt G = tgamma(shape);
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
	register float_mt term;
	
	for(int t = D; t<T-D; t++)
	{
		for(int j = 0; j<N; j++)
		{
			for(int d = 0; d<D; d++)
			{
				for(int i=0; i<N; i++)
				{
					term = alpha[(t-d-1)*N+i]*A[i*N+j]*(dur_dists[j])[d]*prod(&c[t-d],d-1);
					term *= prod(&(ydensity[j])[t-d],d);
					alpha[t*N+j] += term;
				}
			}
		}
		c[t] = 1/sum(&alpha[T*N], N);
		for(int i=0; i<N; i++)
		{
			alpha[t*N + i] *= c[t];
		}
	}
}
static inline void compute_betas(void)
{
	register float_mt term;
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
					term = A[i*N+j]*(dur_dists[j])[d]*prod(&c[t+1], d-1)*beta[(t+d-1)*N+j];
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
	float_mt num;
	float_mt term;
	float_mt den;
	float_mt *temp = malloc(sizeof(float_mt)*T);
	for(int i=0;i<N;i++)
	{
		for(int j=0; j<N; j++)
		{
			num = 0;
			for(int t = D; t<T-D; t++)
			{
				for(int d=0; d<D; d++)
				{
					term = A[i*N+j]*(dur_dists[j])[d]*beta[(t+d)*N+j]*prod(&c[t+1],d-1);
					term *= alpha[t*N+i]*prod(&(ydensity[j])[t+1],d);
					num += term;
				}
			temp[t]=alpha[t*N+i]*beta[t*N+i]/c[t];
			}
			den = sum(&temp[D], T-2*D);
			Ahat[i*N+j] = num/den;
		}
	}
	free(temp);
}
static inline void compute_mu_sigma(void)
{
	float_mt nummu = 0;
	float_mt numsig = 0;
	float_mt den = 0;
	float_mt term;
	float_mt somme;
	for(int j=0; j<N; j++)
	{
		for(int t=D; t<T-D; t++)
		{
			for(int d=0;d<D;d++)
			{
				somme = 0;
				for(int k=0; k<d; k++)
				{
					somme += (y[t-d]-mu(j))*(y[t-d]-mu(j));
				}
				for(int i =0; i<N; i++)
				{
					term = alpha[(t-d+1)*N+i]*A[i*N+j]*(dur_dists[j])[d]*prod(&c[t-d],d-1);
					term *= prod((ydensity[j])[t-d], d-1);
					nummu += term*sum(&y[t-d],d);
					numsig += term*somme;
					den = den+term;
				}
			}
		}
		muhat(j) = nummu/den;
		sigmahat(j) = numsig/den;
	}
}
static inline void compute_duration_params(void)
{

}
static inline void compute_duration_dist_params(void)
{
}
































