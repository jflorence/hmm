#include <stdlib.h>
#include "train.h"
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "aux.h"
#include <math.h>
#include <gsl/gsl_sf_psi.h>
#include "debug.h"


#define MAXCOUNTER 1

#define WALPHA

static int N;
static long long T;
static long D;
static long size_dist;
static long size_dur_dist;
static int scalingstride = 1;
static float_mt step = 0.001;
static float_mt dstep = 1.0;
static float_mt *alpha;
static float_mt *beta;
static float_mt **ydensity;
static float_mt **distribs;
static float_mt **dur_dists;
static delay_mt *y;
static float_mt *A;
static float_mt *Ahat;
static float_mt *mu;
static float_mt *muhat;
static float_mt *sigmahat;
static float_mt *c;
static float_mt *shape;
static float_mt *scale;
static float_mt *dshape;
static float_mt *dscale;
static float_mt *dshapehat;
static float_mt *dscalehat;




static inline void traininit(struct params *p, delays_mt *data, delay_mt ymax);
static void freetrainparams(void);
static inline void initalpha(void);
static inline void initbeta(void);
static inline void compute_distribs(long size_dist, long size_dur_dist);
static inline void compute_gamma_dist(float_mt *dist, int length, float_mt step, float_mt shape, float_mt scale);
static inline void density_of_y();
static inline void compute_alphas(void);
static inline void compute_betas(void);
static inline void compute_A(void);
static inline void compute_mu_sigma(void);
static inline void compute_duration_scale(void);
static inline void compute_duration_shape(void);






void train(struct params *p, delays_mt *data, delay_mt ymax)
{
	float_mt previouslike = -INFINITY;
	float_mt currentlike;
	float_mt epsilon = -1.0; //don't care about this for the moment...
	float_mt Dlike = epsilon+1.0;
	int counter = 0;
	
	traininit(p, data, ymax);	

	while(counter++ < MAXCOUNTER && Dlike>epsilon) //or likelihood...
	{
		dispstr(1,"Starting iteration %d\n", counter);
		compute_distribs(size_dist, size_dur_dist);
		density_of_y();
		initalpha();
		compute_alphas();
		initbeta();
		compute_betas();
		compute_A();
		compute_mu_sigma();
		compute_duration_scale();
		compute_duration_shape();

		//update A, mu, sigma...
		memcpy(A, Ahat, N*N*sizeof(float_mt));
		memcpy(p->A, Ahat, N*N*sizeof(float_mt));
		memcpy(mu, muhat, N*sizeof(float_mt));
		memcpy(p->mu, muhat, N*sizeof(float_mt));
		memcpy(p->sigma, sigmahat, N*sizeof(float_mt));
		//compute gamma params
		for(int i = 0; i<N; i++)
		{
			p->scale[i] = scale[i] = sigmahat[i]/mu[i];
			p->shape[i] = shape[i] = mu[i]*mu[i]/sigmahat[i];
			p->dscale[i] = dscale[i] = dscalehat[i];
			p->dshape[i] = dshape[i] = dshapehat[i];
		}
		//update log likelihood
		currentlike = 0.0;
		for(int i=0; i<T; i++)
		{
			currentlike+=log(c[i]);
		}
		Dlike = fabsl((long double)currentlike-(long double)previouslike);
		previouslike = currentlike;
	}

#ifdef WALPHA
	dispstr(1, "Writing alphas and betas to file...\n");
	write_alpha_beta(alpha, beta, N, T);
#endif


	FILE *file = fopen("c.txt","w");
	for(int i=0; i<T; i++)
	{
		fprintf(file, "%Le\n", c[i]);
	}
	fclose(file);




	freetrainparams();
	
	dispstr(1,"Training complete.\n");
	dispstr(2, "Dlike: %e\n", Dlike);
	dispstr(2, "log-likelihood: %e\n", currentlike);
	return;
}

void compute_distribs(long size_dist, long size_dur_dist)
{
	dispstr(2,"Computing gamma distributions...\n");
	for(int i=0; i<N; i++)
	{
		compute_gamma_dist(distribs[i], size_dist, step, shape[i], scale[i]);
		compute_gamma_dist(dur_dists[i], size_dur_dist, dstep, dshape[i], dscale[i]);	
		//dispstr(3, "i: %d\n", i);	
		//dispstr(4, "size_dist: %ld\n", size_dist);
		//dispstr(4, "step: %Le\n", step);
		//dispstr(4, "shape: %Le\n", shape[i]);
		//dispstr(4, "scale: %Le\n", scale[i]);
		//plot1D(distribs[i], size_dist);
	}
}




//alpha and beta should be stored so that the values for the same t are contiguous.
static inline void initalpha(void)
{
	dispstr(2,"(re)initializing alpha...\n");
	memset(alpha, 0, N*T*sizeof(float_mt));
	alpha[N*(D-1)] = 1.0;
	return;
}



static inline void initbeta(void)
{
	dispstr(2,"(re)initializing beta...\n");
	memset(beta, 0, N*T*sizeof(float_mt));
	for(int i = N*(T-D); i<(T-D+1)*(N);i++)
	{
		beta[i] = 1.0;
	}
}


static inline void compute_gamma_dist(float_mt *dist, int length, float_mt step, float_mt shape, float_mt scale)
{
	float_mt x = 0.0;
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
	dispstr(2,"Computing densities for input vector...\n");
	//we need the proba density for each sample of the y vector, according to the current distributions
	int index;
	for(int j=0; j<N; j++)
	{
		for(int k = 0; k<T; k++)
		{
			index = round(y[k]/step)+1; //this is an index on the vector of the pdf
			(ydensity[j])[k] = (distribs[j])[index];
		}
		//plot1D(distribs[j], size_dist);
	}
	
}


static inline void compute_alphas(void)
{
	dispstr(2,"Computing alphas...\n");
	//Is this the best order for the loops ?
	float_mt term;
	float_mt prodc;
	float_mt prody;
	float_mt product;	
	for(long long t = D; t<T-D; t++)
	{
		for(int j = 0; j<N; j++)
		{
			prodc = 1.0/c[t];
			prody = 1.0;
			for(int d = 0; d<D; d++)
			{
				prodc *= c[t-d];
				prody *= (ydensity[j])[t-d];
				product = prodc*prody*(dur_dists[j])[d+1];
				for(int i=0; i<N; i++)
				{
					term = alpha[(t-d-1)*N+i]*A[i*N+j]*product;
					alpha[t*N+j] += term;
				}
			}
		}
		if(t%scalingstride == 0)
		{
			//c[t] = 1/sum(&alpha[t*N], N);
		}
		for(int j=0; j<N; j++)
		{
			alpha[t*N + j] *= c[t];
		}
	}
}
static inline void compute_betas(void)
{
	dispstr(2,"Computing betas...\n");
	float_mt term;
	float_mt prodc;
	float_mt prody[N];
	float_mt product;
	long index;
	for(long long t = T-D-1; t>=D; t--)
	{
		for(int i =0; i<N; i++)
		{
			index = t*N+i;
			prodc = 1.0;
			for(int j=0; j<N; j++)
			{
				prody[j] = 1.0;
			}
			for(int d=0; d<D; d++)
			{
				if(d<D-1)
					prodc *= c[t+1+d];
				for(int j = 0; j<N; j++)
				{
					prody[j] *= (ydensity[j])[t+1+d];
					product = prodc*prody[j]*(dur_dists[j])[d+1]*beta[(t+d+1)*N+j];
					term = A[i*N+j]*product;
					beta[index] += term;
				}
			}
			beta[index] *= c[t];
		}
	}
}
static inline void compute_A(void)
{
	dispstr(2,"Computing A...\n");
	float_mt num;
	float_mt term;
	float_mt den;
	float_mt prodc;
	float_mt prody;
	float_mt product;
	float_mt somme;
	float_mt *temp = malloc(sizeof(float_mt)*T);
	if(temp == NULL)
	{
		fprintf(stderr, "Couldn't malloc\n");
		exit(EXIT_FAILURE);
	}
	for(int i=0;i<N;i++)
	{
		for(int j=0; j<N; j++)
		{
			num = 0.0;
			for(long long t = D; t<T-D; t++)
			{
				prodc = 1.0;
				prody = 1.0;
				somme = 0.0;
				for(int d=0; d<D; d++)
				{
					if(d<D-1)
						prodc *= c[t+1+d];
					prody *= (ydensity[j])[t+d+1];
					product = prody*prodc*(dur_dists[j])[d+1]*beta[(t+d+1)*N+j]*A[i*N+j];
					somme += product;
					term = alpha[t*N+i]*product;
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
	dispstr(2,"Computing mu and sigma...\n");
	float_mt nummu = 0.0;
	float_mt numsig = 0.0;
	float_mt den = 0.0;
	float_mt term;
	float_mt somme;
	float_mt somme2;
	float_mt prodc;
	float_mt prody;
	float_mt product;
	for(int j=0; j<N; j++)
	{
		for(long t=D; t<T-D; t++)
		{
			prodc = 1.0/c[t];
			prody = 1.0;
			somme = 0;
			somme2 = 0;
			for(int d=0;d<D;d++)
			{
				somme += (y[t-d]-mu[j])*(y[t-d]-mu[j]);
				somme2 += y[t-d];
				prodc *= c[t-d];
				prody *= (ydensity[j])[t-d];
				
				product = prodc*prody*(dur_dists[j])[d+1];
				for(int i =0; i<N; i++)
				{
					term = alpha[(t-d-1)*N+i]*A[i*N+j]*product;
					nummu += term*somme2;
					numsig += term*somme;
					den = den+term;
				}
			}
		}
		muhat[j] = nummu/den;
		sigmahat[j] = numsig/den;
	}
}
static inline void compute_duration_scale(void)
{
	dispstr(2,"Computing duration scale param...\n");
	float_mt num;
	float_mt term;
	float_mt somme;
	float_mt den;
	float_mt prodc;
	float_mt prody;
	float_mt product;
	for(int j=0;j<N;j++)
	{
		num = 0.0;
		for(long long t=D;t<T-D;t++)
		{
			prodc = 1.0/c[t];
			prody = 1.0;
			for(int d=0;d<D;d++)
			{
				prodc *= c[t-d];
				prody *= (ydensity[j])[t-d];
				product = prodc*prody*(dur_dists[j])[d+1]*(d+1)*beta[t*N+j];
				for(int i=0; i<N; i++)
				{
					term = alpha[(t-d-1)*N+i]*A[i*N+j]*product;
					num += term;
				}
			}
		}
		somme = 0.0;
		for(long long k=D; k<T-D; k++)
		{
			somme += alpha[k*N+j]*beta[k*N+j]/c[k];
		}
		den = dshape[j]*somme;
		dscalehat[j] = num/den;
	}
}
static inline void compute_duration_shape(void)
{
	dispstr(2,"Computing duration shape param...\n");
	float_mt C;
	float_mt digamma;
	float_mt polygamma1;
	float_mt num;
	float_mt term;
	float_mt somme;
	float_mt *rate = malloc(sizeof(float_mt)*N);
	
	float_mt prodc;
	float_mt prody;
	float_mt product;

	if(rate == NULL)
	{
		fprintf(stderr, "Couldn't malloc\n");
		exit(EXIT_FAILURE);
	}
	for(int i=0; i<N; i++)
	{
		rate[i] = 1.0/dscale[i];
	}
	for(int j=0; j<N; j++)
	{
		digamma = gsl_sf_psi_n(0, (double) dshape[j]);
		polygamma1 = gsl_sf_psi_n(1, (double) dshape[j]);
		num = 0.0;
		for(long long t=D; t<T-D; t++)
		{
			prodc = 1.0/c[t];
			prody = 1.0;
			for(int d=0; d<D; d++)
			{
				prodc *= c[t-d];
				prody *= (ydensity[j])[t-d];
				product = prodc*prody*(dur_dists[j])[d+1]*log(rate[j]*(d+1))*beta[t*N+j];
				for(int i=0; i<N; i++)
				{
					term = product*alpha[(t-d-1)*N+i]*A[i*N+j];
					num += term;
				}
			}
		}
		somme = 0.0;
		for(long long k=D; k<T-D; k++)
		{
			somme += alpha[k*N+j]*beta[k*N+j]/c[k];
		}
		C = num/somme;
		dshapehat[j] = dshape[j]-(digamma-C)/polygamma1;
	
		if(dshapehat[j] <= 0.0)
		{
			dshapehat[j] = 0.01;
		}
	}
}




static void freetrainparams(void)
{
	free(alpha);
	free(beta);
	for(int i = 0; i<N; i++)
	{
		free(distribs[i]);
		free(dur_dists[i]);
		free(ydensity[i]);
	}
	free(distribs);
	free(dur_dists);
	free(ydensity);
	free(A);
	free(c);
	free(mu);
	free(muhat);
	free(sigmahat);
	free(shape);
	free(scale);
	free(dshape);
	free(dscale);
	free(dshapehat);
	free(dscalehat);


}



static inline void traininit(struct params *p, delays_mt *data, delay_mt ymax)
{
	dispstr(1,"Initializing training algo...\n");
	N = p->N;
	T = data->length;
	D = p->D;
	y = data->delay;

	alpha = malloc(N*T*sizeof(float_mt));
	beta = malloc(N*T*sizeof(float_mt));
	if(alpha == NULL || beta == NULL)
	{
		fprintf(stderr, "Couldn't malloc\n");
		exit(EXIT_FAILURE);
	}




	//T = 2*D+10;

	printf("t: %lld\n", T);
	printf("d: %ld", D);



	assert(T>2*D);
	initalpha();
	initbeta();

	size_dist = ceil(ymax/step)+1;
	size_dur_dist = D+1;
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
	c = malloc(T*sizeof(float_mt));
	Ahat = malloc(N*N*sizeof(float_mt));
	mu = malloc(sizeof(float_mt)*N);
	muhat = malloc(sizeof(float_mt)*N);
	sigmahat = malloc(sizeof(float_mt)*N);
	shape = malloc(sizeof(float_mt)*N);
	scale = malloc(sizeof(float_mt)*N);
	dshape = malloc(sizeof(float_mt)*N);
	dscale = malloc(sizeof(float_mt)*N); 
	dshapehat = malloc(sizeof(float_mt)*N); 
	dscalehat = malloc(sizeof(float_mt)*N); 
	if(A == NULL || Ahat == NULL || mu == NULL || muhat == NULL || sigmahat == NULL ||
		dshape == NULL || dscale == NULL || dshapehat == NULL || dscalehat == NULL || 
		shape == NULL || scale == NULL)
	{
		fprintf(stderr, "Couldn't malloc\n");
		exit(EXIT_FAILURE);
	}
	memcpy(A, p->A, N*N*sizeof(float_mt));
	for(long long i = 0; i<T; i++)
		c[i] = 1.0;
	for(int i = 0; i<N;i++)
	{
		mu[i] = p->mu[i];
		scale[i] = p->scale[i];
		shape[i] = p->shape[i];
		dscale[i] = p->dscale[i];
		dshape[i] = p->dshape[i];
	}	
	return;	
}




